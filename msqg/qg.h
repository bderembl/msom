/**
   multiple scale quasi-geostrophic model 
*/
#define MODE_PV_INVERT 0

double * dhf;
double * dhc;
double * idh0;
double * idh1;

#include "layer.h"
#include "predictor-corrector.h"
#include "poisson.h"
#include "poisson_layer.h"
#include "timestep.h"
#include "bcg.h"

// for mkdir
#include <sys/stat.h>
#include <sys/types.h>

// layered variables
scalar * qol  = NULL; // vorticity on layers
scalar * pol  = NULL; // stream function on layers
scalar * zetal = NULL; // relative vorticity
scalar * qofl = NULL; // filter mean
//scalar * qosl = NULL; // filter mean

scalar * qom  = NULL; // vorticity on vertical modes
scalar * pom  = NULL; // stream function on modes

// large scale variables (inversion matrices, etc)
scalar * ppl  = NULL; // large scale stream function
scalar * strl = NULL; // vertical stretching (f^2/N^2)
// If PV inversion is done with the tri-diagonal method, 
// then no need to refine cl2m,cm2l,iBul, Frl, Ro
scalar * cl2m  = NULL;
scalar * cm2l  = NULL;
scalar * iBul = NULL; // inverse burger number
scalar * Frl = NULL;
scalar * ptracersl = NULL;   // ptracer field
scalar * ptr_relaxl = NULL;  // ptracer relaxation field
scalar * ptr_srcl = NULL;   // ptracer source term


scalar Ro[];
scalar Rd[];
scalar sig_filt[];
scalar sig_lev[];
scalar topo[];
int flag_topo;
int lsmin, lsmax; 
double afilt = 10.;  // filter size = afilt*Rd
double Lfmax = 1.e10;  // max filter length scale

// temporary psi field with samem BC as p
scalar * tmpl = NULL;

scalar * evolving = NULL;

mgstats mgpsi;

int nl = 1;

double Re = 0.; // reynolds number
double Re4 = 0.; // bihormonic reynolds number
double iRe = 0.0; // inverse reynolds number
double iRe4 = 0.0; // inverse bihormonic reynolds number
double Ekb = 0.0; // Ekman number (bottom)
double Eks = 0.0; // Ekman number (surface)
double Rom = 0.0; // Mean Rossby number
double Frm[1000]; // Mean Froude number
double dhu[1000]; // user input dh
double upg[1000] = {0};  // background U
double vpg[1000] = {0};  // background V
double beta = 0.5;
double tau0 = 0.; // wind stress curl
double sbc = 0.;  // doubly periodic: -1, free slip: 0: (OK), so slip: big number (TO BE FINISHED)

double tend = 1; // end time
double dtflt = -1; // Delat T filtering
double dtout = 1; // Delat T output
double ptr_r[1000] = {0};  // inverse relaxation time scale for passive tracers

int nbar = 0;
int ediag = -1;  // ediag = -1: no ediag, 0: psi*dqdt, 1: (psi+pg)*dqdt
int varRo = 0;   // varRo = 1: variable Rossby number (multiple scale mode)
int l_tmp = 0;   // global layer variable for periodic BC (to be replaced by _layer in the new layer framework)
int nptr = 0;   // number of passive tracers

char dpath[80]; // name of output dir

#if MODE_PV_INVERT
#include "eigmode.h"
#endif

trace
void invertq(scalar * pol, scalar * qol)
{
#if MODE_PV_INVERT
  // convert to modes: matrix multiplication cl2m*qo
  foreach(){
    // outer loop on modes
    for (int m = 0; m < nl ; m++) {
      scalar qm  = qom[m];
      qm[] = 0.;
      
      // inner loop on layers
      for (int l = 0; l < nl ; l++) {
        scalar qo  = qol[l];
        scalar l2m = cl2m[m*nl+l];
        qm[] += l2m[]*qo[];
        }      
      }
    }

  boundary(qom);
//  boundary(pom);
  // elliptic solver on modal pressures
  for (int l = 0; l < nl ; l++) {
    scalar qm  = qom[l];
    scalar pm  = pom[l];
    scalar iBu  = iBul[l];
    mgpsi = poisson (pm, qm, lambda = iBu, tolerance = 1e-3);
  }

  // back to layers : matrix multiplication cm2l*pm
  foreach(){
    // outer loop on layers
    for (int l = 0; l < nl ; l++) {
      scalar po  = pol[l];
      po[] = 0.;
      
      // inner loop on modes
      for (int m = 0; m < nl ; m++) {
        scalar pm  = pom[m];
        scalar m2l = cm2l[l*nl+m];
        po[] += m2l[]*pm[];
      }
    }
  }
#else
  mgpsi = poisson_layer (pol, qol, strl = strl, tolerance = 1e-3);
#endif

  boundary(pol);
}

/**
   relative vorticity
*/

#define laplacian(po) (po[1] + po[-1] + po[0,1] + po[0,-1] - 4*po[])/(sq(Delta))

trace
void comp_del2(scalar * pol, scalar * zetal, double add, double fac)
{
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar po = pol[l];
      scalar zeta = zetal[l];
      zeta[] = add*zeta[] + fac*laplacian(po);
    }

  boundary(zetal);
  
  /**
   Slip BC*/
  if (sbc > 0) {
    for (int l = 0; l < nl ; l++) {
      scalar po = pol[l];
      scalar zeta = zetal[l];
      foreach_boundary(left) {
        zeta[ghost] = sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
      foreach_boundary(right) {
        zeta[ghost] = sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
      foreach_boundary(top) {
        zeta[ghost] = sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
      foreach_boundary(bottom) {
        zeta[ghost] = sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
    }
  }

}

trace
void comp_stretch(scalar * pol, scalar * stretchl, double add, double fac)
{

  foreach() {

    if (nl > 1){
      // upper layer
      int l = 0;
      scalar po_1 = pol[l];
      scalar po_2 = pol[l+1];
      scalar stretch = stretchl[l];
      scalar s1 = strl[l];

      stretch[] = add*stretch[] + fac*s1[]*(po_2[] - po_1[])*idh1[l];

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        scalar po_0 = pol[l-1];
        scalar po_1 = pol[l];
        scalar po_2 = pol[l+1];
        scalar stretch = stretchl[l];
        scalar s0 = strl[l-1];
        scalar s1 = strl[l];

        stretch[] = add*stretch[] + fac*(s0[]*(po_0[] - po_1[])*idh0[l] + s1[]*(po_2[] - po_1[])*idh1[l]);
      }
      // lower layer
      l = nl-1;
      scalar po_0 = pol[l-1];
      po_1 = pol[l];
      stretch = stretchl[l];
      scalar s0 = strl[l-1];

      stretch[] = add*stretch[] + fac*s0[]*(po_0[] - po_1[])*idh0[l];
    }
    else{
      scalar stretch = stretchl[0];
      stretch[] = 0.;
    }
  }

  boundary(stretchl);
}
/**
   Arakawa jacobian: Energy/esntrophy conserving scheme
   This macro computes -J(p,q)
*/

#define jacobian(po,qo) ((( qo[1, 0 ]-qo[-1, 0])*(po[0, 1]-po[ 0 ,-1]) \
     +(qo[0 ,-1]-qo[ 0 ,1])*(po[1, 0]-po[-1, 0 ]) \
     + qo[1, 0 ]*( po[1,1 ] - po[1,-1 ]) \
     - qo[-1, 0]*( po[-1,1] - po[-1,-1]) \
     - qo[ 0 ,1]*( po[1,1 ] - po[-1,1 ]) \
     + qo[0 ,-1]*( po[1,-1] - po[-1,-1]) \
     + po[ 0 ,1]*( qo[1,1 ] - qo[-1,1 ]) \
     - po[0 ,-1]*( qo[1,-1] - qo[-1,-1]) \
     - po[1, 0 ]*( qo[1,1 ] - qo[1,-1 ]) \
     + po[-1, 0]*( qo[-1,1] - qo[-1,-1]))\
    /(12.*Delta*Delta))

/**
   - beta v
   (note the minus sign because it is in the rhs)
*/

#define beta_effect(po) (beta*(po[-1] - po[1])/(2*Delta))

/**
   Compute velocity at faces
*/

trace
void comp_vel(const scalar po, face vector uf)
{
  struct { double x, y; } f = {-1.,1.};
  foreach_face() {
    uf.x[] = f.x*0.25*(po[0,1] - po[0,-1] + po[-1,1] - po[-1,-1])/Delta;
  }
  boundary ((scalar *){uf});
}

// todo: change names qo, zeta!
trace
double advection_pv(scalar * qol, scalar * qotl, scalar * pol, scalar * dqol, double dtmax)
{

  foreach() {
    double ju,jd;

    if (nl > 1){

      // upper layer
      int l = 0;
      scalar qo  = qol[l];
      scalar po  = pol[l];
      scalar pp  = ppl[l];
      scalar qot  = qotl[l];
      scalar dqo = dqol[l];
      scalar po2  = pol[l+1];
      scalar pp2  = ppl[l+1];
      scalar s1 = strl[l];

#if ENERGY_CONSERV
      jd =  jacobian(pp, po2) + jacobian(po, pp2);
      dqo[] += jacobian(po, qot) + jacobian(pp, qo) + beta_effect(po) + s1[]*jd*idh1[l];
#else
      jd = jacobian(po, po2) + jacobian(pp, po2) + jacobian(po, pp2);
      dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s1[]*jd*idh1[l];
#endif

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        qo  = qol[l];
        po  = pol[l];
        qot = qotl[l];
        pp  = ppl[l];
        dqo = dqol[l];
        po2  = pol[l+1];
        pp2  = ppl[l+1];
        scalar s0 = strl[l-1];
        scalar s1 = strl[l];

        ju = -jd;
#if ENERGY_CONSERV
        jd = jacobian(pp, po2) + jacobian(po, pp2);
        dqo[] += jacobian(po, qot) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju*idh0[l] + s1[]*jd*idh1[l];
#else
        jd = jacobian(po, po2) + jacobian(pp, po2) + jacobian(po, pp2);
        dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju*idh0[l] + s1[]*jd*idh1[l];
#endif
      }

      // lower layer
      l = nl-1;

      qo  = qol[l];
      po  = pol[l];
      qot  = qotl[l];
      pp  = ppl[l];
      dqo = dqol[l];
      scalar s0 = strl[l-1];

      ju = -jd;
#if ENERGY_CONSERV
      dqo[] += jacobian(po, qot) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju*idh0[l];
#else
      dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju*idh0[l];
#endif
    }
    else{
      scalar dqo = dqol[0];
      dqo[] = 0.;
    }
  }

  // compute dtmax
   for (int l = 0; l < nl ; l++) {
     face vector uf[];
     scalar po = pol[l];
     comp_vel(po, uf);
     dtmax = timestep (uf, dtmax);
     po = ppl[l];
     comp_vel(po, uf);
     dtmax = timestep (uf, dtmax);
   }
  return dtmax;
}

trace
void comp_q(scalar * pol, scalar * qol)
{
  comp_del2  (pol, qol, 0., 1.);
  comp_stretch(pol, qol, 1., 1.);
  // TODO: not the rigght BC if partial or no slip
  boundary(qol);
}


trace
void dissip  (scalar * zetal, scalar * dqol)
{
  comp_stretch(zetal, dqol, 1., iRe);

  comp_del2(zetal, tmpl, 0., 1.);

  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar dqo = dqol[l];
      scalar p4 = tmpl[l];
      dqo[] += p4[]*iRe;
    }

  comp_stretch(tmpl, dqol, 1., iRe4);
  comp_del2(tmpl, dqol, 1., iRe4);
}

/**
   Bottom and surface Ekman Friction
*/

trace
void ekman_friction  (scalar * zetal, scalar * dqol)
{
  foreach() {
    scalar dqos = dqol[0];
    scalar zetas = zetal[0];

    scalar dqob = dqol[nl-1];
    scalar zetab = zetal[nl-1];
    dqos[] -= Eks/(Rom*2*dhf[0])*zetas[];
    dqob[] -= Ekb/(Rom*2*dhf[nl-1])*zetab[];
  }
}

/**
   surface forcing
*/

trace
void surface_forcing  (scalar * dqol)
{
  scalar dqo = dqol[0];
  foreach()
    dqo[] -= tau0/(Rom*dhf[0])*sin(2*pi*y/L0)*sin(pi*y/L0);
}

/**
   Bottom topography
*/

trace
void bottom_topography  (scalar * pol, scalar * dqol)
{
  foreach() {
    scalar dqo = dqol[nl-1];
    scalar po  = pol[nl-1];
    dqo[] += jacobian(po, topo)/(Rom*dhf[nl-1]);
  }
}


trace
void time_filter (scalar * qol, scalar * qo_mel, double dt)
{

  double tau_f = 20;
  double alpha_f = dt/tau_f;


  scalar qo_me[];
  scalar qo[];
  foreach(){
    for (qo_me, qo in qo_mel, qol){
      qo_me[] = (1-alpha_f)*qo_me[] + alpha_f*qo[];
//      qo[] = qo[] - qo_me[];
    }
  }
}

trace
void wavelet_filter(scalar *qol, scalar * pol, scalar * qofl, double dtflt, int nbar)
{

  /* save q values */
//  scalar * qosl = list_clone(qol);
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar qo  = qol[l];
      scalar tmp  = tmpl[l];
      tmp[] = qo[];
    }

  invertq(pol,qol);
  
  for (scalar po in pol) {
    scalar w[];
    w[top] = 0;
    w[bottom] = 0;
    w[right] = 0;
    w[left] = 0;


    wavelet(po,w);
    for (int l = 0; l <= depth(); l++) {
      foreach_level (l)
        w[] *= sig_lev[];
      boundary_level ({w}, l);
    }
    inverse_wavelet (po, w);

  }
  
  comp_q(pol,qol);
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar qo  = qol[l];
      scalar qof  = qofl[l];
      scalar tmp  = tmpl[l];

      qof[] = (qof[]*nbar + (tmp[] - qo[])/dtflt)/(nbar+1);
      //qof[] = (tmp[] - qo[]);
    }
  // for energy diag: restore qo to prefiltered value
  if (dtflt < 0.0){
    list_copy_deep (tmpl, qol, nl);
  }

  boundary(qofl);
  nbar++;
  
}

/**
   ## time stepping routines
   We use the predictor corrector implementation */

static void advance_qg (scalar * output, scalar * input,
                        scalar * updates, double dt)
{
  foreach() {
    for (int l = 0; l < nl ; l++) {
      scalar qi = input[l];
      scalar qo = output[l];
      scalar dq = updates[l];
      qo[] = qi[] + dq[]*dt;
    }
  }
  boundary(output);
}

double update_qg (scalar * evolving, scalar * updates, double dtmax)
{
  foreach()
    for (scalar s in updates)
      s[] = 0.;

  invertq(pol, evolving);
  comp_del2(pol, zetal, 0., 1.0);
  dtmax = advection_pv(zetal, qol, pol, updates, dtmax);
  dissip(zetal, updates);
  ekman_friction(zetal, updates);
  surface_forcing(updates);
  if (flag_topo)
    bottom_topography(pol,updates);
  
  return dtmax;
}

/**
   Filter
*/
event filter (t = dtflt; t <= tend+1e-10;  t += dtflt) {
  fprintf(stdout,"Filter solution\n");
  wavelet_filter ( qol, pol, qofl, dtflt, nbar)
}
/**
   Ptracers
*/

event tracer_advection (i++,last) {

  if (nptr>0) 
    for (int l = 0; l < nl ; l++) {
      face vector uf[];
      scalar po = pol[l];
      comp_vel(po, uf);
      
      scalar * list_ptr_lev = NULL;
      scalar * list_src_lev = NULL;
      for (int nt = 0; nt < nptr ; nt++) {
        scalar ptracers = ptracersl[l*nptr + nt];
        scalar ptr_relax = ptr_relaxl[l*nptr + nt];
        scalar ptr_src = ptr_srcl[l*nptr + nt];
        list_ptr_lev = list_append (list_ptr_lev, ptracers);
        foreach()
          ptr_src[] = ptr_r[nt]*(ptr_relax[] - ptracers[]);
        
        list_src_lev = list_append (list_src_lev, ptr_src);
      }
      boundary(list_src_lev);
      advection (list_ptr_lev, uf, dt, list_src_lev);
      free(list_ptr_lev);
      free(list_src_lev);
    }
}


/**********************************************************************
*                       End of dynamical core                         *
***********************************************************************/

/**
   Read input parameters
 */

void trim_whitespace(char* s) {
  const char* d = s;
  do {
    while (*d == ' ')
      ++d;
  } while (*s++ = *d++);
}


void str2array(char *tmps2, double *array){
  char* p;
  trim_whitespace(tmps2);
  int len = strlen(tmps2);
  char tmps3[len];
  strcpy(tmps3, tmps2); //needed in case there is an empty line in params.in
  int n = 0;
  p = strtok(tmps3,"[,]");
  while (p != NULL){
    array[n] = atof(p);
    p = strtok(NULL, ",");
    n += 1;
  }
}

void read_params(char* path2file)
{
  FILE * fp;
  if ((fp = fopen(path2file, "rt"))) {
    char tempbuff[300];
    while(fgets(tempbuff,300,fp)) {
      char* tmps1 = strtok(tempbuff, "=");
      char* tmps2 = strtok(NULL, "=");
      trim_whitespace(tmps1);
      if      (strcmp(tmps1,"N")    ==0) { N     = atoi(tmps2); }
      else if (strcmp(tmps1,"nl")   ==0) { nl    = atoi(tmps2); }
      else if (strcmp(tmps1,"ediag")==0) { ediag = atoi(tmps2); }
      else if (strcmp(tmps1,"varRo")==0) { varRo = atoi(tmps2); }
      else if (strcmp(tmps1,"nptr") ==0) { nptr  = atoi(tmps2); }
      else if (strcmp(tmps1,"L0")   ==0) { L0    = atof(tmps2); }
      else if (strcmp(tmps1,"Rom")  ==0) { Rom   = atof(tmps2); }
      else if (strcmp(tmps1,"Ekb")  ==0) { Ekb   = atof(tmps2); }
      else if (strcmp(tmps1,"Eks")  ==0) { Eks   = atof(tmps2); }
      else if (strcmp(tmps1,"tau0") ==0) { tau0  = atof(tmps2); }
      else if (strcmp(tmps1,"Re")   ==0) { Re    = atof(tmps2); }
      else if (strcmp(tmps1,"Re4")  ==0) { Re4   = atof(tmps2); }
      else if (strcmp(tmps1,"sbc")  ==0) { sbc   = atof(tmps2); }
      else if (strcmp(tmps1,"beta") ==0) { beta  = atof(tmps2); }
      else if (strcmp(tmps1,"afilt")==0) { afilt = atof(tmps2); }
      else if (strcmp(tmps1,"Lfmax")==0) { Lfmax = atof(tmps2); }
      else if (strcmp(tmps1,"DT")   ==0) { DT    = atof(tmps2); }
      else if (strcmp(tmps1,"tend") ==0) { tend  = atof(tmps2); }
      else if (strcmp(tmps1,"dtout")==0) { dtout = atof(tmps2); }
      else if (strcmp(tmps1,"dtflt")==0) { dtflt = atof(tmps2); }
      else if (strcmp(tmps1,"CFL")  ==0) { CFL   = atof(tmps2); }
      else if (strcmp(tmps1,"Fr")   ==0) { str2array(tmps2, Frm);}
      else if (strcmp(tmps1,"dh")   ==0) { str2array(tmps2, dhu);}
      else if (strcmp(tmps1,"upg")  ==0) { str2array(tmps2, upg);}
      else if (strcmp(tmps1,"vpg")  ==0) { str2array(tmps2, vpg);}
      else if (strcmp(tmps1,"ptr_r")==0) { str2array(tmps2, ptr_r);}
//      printf("%s => %s\n", tmps1, tmps2);
    }
    fclose(fp);
  } else {
    fprintf(stdout, "file %s not found\n", path2file);
    exit(0);
  }
  if (Re  == 0) iRe  = 0.; else iRe  =  1/Re;
  if (Re4 == 0) iRe4 = 0.; else iRe4 = -1/Re4;

  /**
     Viscosity CFL = 0.5
   */
  if (Re  != 0) DT = 0.5*min(DT,sq(L0/N)*Re/4.);
  if (Re4 != 0) DT = 0.5*min(DT,sq(sq(L0/N))*Re4/32.);

  fprintf(stdout, "Config: N = %d, nl = %d, L0 = %g\n", N, nl, L0);
}

/**
   Create output directory and copy input parameter file for backup
*/
void create_outdir()
{
  if (pid() == 0) {
    for (int i=1; i<10000; i++) {
      sprintf(dpath, "outdir_%04d/", i);
      if (mkdir(dpath, 0777) == 0) {
        fprintf(stdout,"Writing output in %s\n",dpath);
        break;
      }
    }
  }
@if _MPI
  MPI_Bcast(&dpath, 80, MPI_CHAR, 0, MPI_COMM_WORLD);
@endif
}

void backup_config()
{
  fprintf(stdout, "Backup config\n");
  char ch;
  char name[80];
  sprintf (name,"%sparams.in", dpath);
  FILE * source = fopen("params.in", "r");
  FILE * target = fopen(name, "w");
  while ((ch = fgetc(source)) != EOF)
    fputc(ch, target);
  fclose(source);
  fclose(target);

  sprintf (name,"%ssig_filt.bas", dpath);
  FILE * fp = fopen (name, "w");
  output_matrixl ({sig_filt}, fp);
  fclose(fp);

#if MODE_PV_INVERT
  sprintf (name,"%siBu.bas", dpath);
  fp = fopen (name, "w");
  output_matrixl (iBul, fp);
  fclose(fp);
#else
  sprintf (name,"%srdpg_%dl_N%d.bas", dpath, nl,N);
  fp = fopen (name, "w");
  output_matrixl ({Rd}, fp);
  fclose(fp);
#endif

  sprintf (name,"%spsipg_%dl_N%d.bas", dpath, nl,N);
  fp = fopen (name, "w");
  output_matrixl (ppl, fp);
  fclose(fp);

  sprintf (name,"%sfrpg_%dl_N%d.bas", dpath, nl,N);
  fp = fopen (name, "w");
  output_matrixl (Frl, fp);
  fclose(fp);

  float dh[nl]; // float instead of double
  for (int l = 0; l < nl ; l++)
    dh[l] = dhf[l];

  sprintf (name,"%sdh_%dl.bin", dpath, nl);
  fp = fopen (name, "w");
  fwrite(&dh, sizeof(float), nl, fp);
  fclose(fp);
}

void set_vars()
{

  int bc_type = 0;

  if (sbc == -1){
    periodic(right);
    periodic(top);
    bc_type = -2; // any negative number < -1
  }

  fprintf(stdout,"Create main variables .. ");
  pol   = create_layer_var(pol,nl,bc_type);
  qol   = create_layer_var(qol,nl,bc_type);
  ppl   = create_layer_var(ppl,nl,bc_type);
  zetal = create_layer_var(zetal,nl,bc_type);
  tmpl  = create_layer_var(tmpl,nl,bc_type);
  Frl   = create_layer_var(Frl,nl,bc_type+1); // periodic or neumann
  strl = create_layer_var(strl,nl,bc_type+1); // periodic or neumann
  qofl  = create_layer_var(qofl,nl,bc_type);
#if MODE_PV_INVERT
  pom = create_layer_var(pom,nl,bc_type);
  qom = create_layer_var(qom,nl,bc_type);
  iBul  = create_layer_var(iBul,nl,bc_type+1);
#endif

  if (nptr > 0){
    ptracersl = create_layer_var(ptracersl,nl*nptr,bc_type+1); // periodic or neumann
    ptr_srcl = create_layer_var(ptr_srcl,nl*nptr,bc_type+1); // periodic or neumann
    ptr_relaxl = create_layer_var(ptr_relaxl,nl*nptr,bc_type+1); // periodic or neumann
  }

  /**
     Mode to layer inversion matrices (dimesnion: $nl^2$)
  */
#if MODE_PV_INVERT
  int nl2 = nl*nl;
  cl2m = create_layer_var(cl2m,nl2,bc_type+1);
  cm2l = create_layer_var(cm2l,nl2,bc_type+1);
#endif
  evolving = qol;
    
  /**
     Default variables:
     Layer thicknesses dhf (dhc is computed after)
     Froude number Fr
     Rossby number Ro
  */
  dhc = malloc ((nl-1)*sizeof(double));
  dhf = malloc (nl*sizeof(double));
  idh0 = malloc (nl*sizeof(double));
  idh1 = malloc (nl*sizeof(double));
  for (int l = 0; l < nl; l++)
    dhf[l] = dhu[l];

  foreach()
    for (int l = 0; l < nl-1 ; l++) {
      scalar Fr = Frl[l];
      Fr[] =  Frm[l];
    }

  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar pp = ppl[l];
      pp[] =  vpg[l]*x - upg[l]*y;
    }

  foreach(){
    Ro[] = Rom; 
    Rd[] = 1.; 
    topo[] = 0.;
  }

  /**
     We overload the default 'advance' and 'update' functions of the
     predictor-corrector scheme and (TODO) setup the prolongation and
     restriction methods on trees. */

  advance = advance_qg;
  update = update_qg;
  fprintf(stdout,"ok\n");
}

event defaults (i = 0){
  set_vars();
}

void set_const() {

/**
   Layer thickness and large scale variables
*/
  fprintf(stdout, "Read input files:\n");

  char name[80];
  FILE * fp;
  sprintf (name,"dh_%dl.bin", nl);
  if ((fp = fopen (name, "r"))) {
    float dh[nl];
    fread(&dh, sizeof(float), nl, fp);
    fclose(fp);
    for (int l = 0; l < nl ; l++)
      dhf[l] = dh[l];
    fprintf(stdout, "%s .. ok\n", name);
  }

  sprintf (name,"psipg_%dl_N%d.bas", nl,N);
  if ((fp = fopen (name, "r"))) {
    input_matrixl (ppl, fp);
    fclose(fp);
    fprintf(stdout, "%s .. ok\n", name);
  }

  sprintf (name,"frpg_%dl_N%d.bas", nl,N);
  if ((fp = fopen (name, "r"))) {
    input_matrixl (Frl, fp);
    fclose(fp);
    fprintf(stdout, "%s .. ok\n", name);
  }

  sprintf (name,"rdpg_%dl_N%d.bas", nl,N);
  if ((fp = fopen (name, "r"))) {
    input_matrixl ({Rd}, fp);
    fclose(fp);
    fprintf(stdout, "%s .. ok\n", name);
  }

  sprintf (name,"topo.bas", nl,N);
  if ((fp = fopen (name, "r"))) {
    flag_topo = 1;
    input_matrixl ({topo}, fp);
    fclose(fp);
    fprintf(stdout, "%s .. ok\n", name);
  }

  /**
     Sanity checks
   */
  for (int l = 0; l < nl ; l++) {
    if (dhf[l] == 0){
      fprintf(stdout, "thickness = 0: aborting\n");
      fprintf(stdout, "Check the definition of dh in params.in\n");
      exit(0);
      }
  }
  foreach()
    for (int l = 0; l < nl-1 ; l++) {
      scalar Fr = Frl[l];
      if (Fr[] == 0){
        fprintf(stdout, "Fr = 0: aborting\n");
        fprintf(stdout, "Check the definition of Fr in params.in\n");
        exit(0);
      }
    }
  if (Rom <= 0){
    fprintf(stdout, "Rom <= 0: aborting\n");
    exit(0);
  }

  /**
     Compute layer metrics
   */
  for (int l = 0; l < nl-1; l++)
    dhc[l] = 0.5*(dhf[l] + dhf[l+1]);

  idh0[0] = 0.;
  idh1[0] = 1./(dhc[0]*dhf[0]);
  for (int l = 1; l < nl-1 ; l++) {
    idh0[l] = 1./(dhc[l-1]*dhf[l]);
    idh1[l] = 1./(dhc[l]*dhf[l]);
  }
  idh0[nl-1] = 1./(dhc[nl-2]*dhf[nl-1]);
  idh1[nl-1] = 0.;

  /**
     Adjust variable rossby number
   */
  if (varRo > 0)
    foreach()
      Ro[] = Rom/(1 + Rom*beta*(y-0.5*L0));
  else
    foreach()
      Ro[] = Rom;

  /**
     Compute vertical stretching coef matrix (inverse burger number squared)
   */

  foreach()
    for (int l = 0; l < nl-1 ; l++) {
      scalar Fr = Frl[l];
      scalar s = strl[l];      
      s[] = sq(Fr[]/Ro[]);
    }

  /**
     compute filter length scale and wavelet coeffs*/
#if MODE_PV_INVERT
  eigmod(dhf, dhc, Ro, Frl, cl2m, cm2l, iBul);

  scalar iRd = iBul[1];
  foreach()
    sig_filt[] = min(afilt*sqrt(-1/iRd[]),Lfmax);
#else
  foreach()
    sig_filt[] = min(afilt*Rd[],Lfmax);
#endif

  restriction({sig_filt});

  // low pass filter
  for (int l = depth(); l >= 0; l--) {
    foreach_level (l) {
      double ref_flag = 0;
      if (l < depth())
        foreach_child()
          ref_flag += sig_lev[];
      if (ref_flag > 0)
        sig_lev[] = 1;
      else
        if (sig_filt[] > 2*Delta)
          sig_lev[] = 0;
        else if (sig_filt[] <= 2*Delta && sig_filt[] > Delta)
          sig_lev[] = 1-(sig_filt[]-Delta)/Delta;
        else
          sig_lev[] = 1;
    }
    boundary_level ({sig_lev}, l);
}

  // high pass filter
  for (int l = depth(); l >= 0; l--) {
    foreach_level (l)
      sig_lev[] = 1 - sig_lev[];
    boundary_level ({sig_lev}, l);  
  }

  comp_q(pol,qol);

  if (nptr>0)
    list_copy_deep (ptracersl, ptr_relaxl, nl*nptr);

  /**
     BC for all fields. If periodic BC, we also adjust the large-scale stream
     function that is not periodic
  */
  boundary (all);

  if (sbc == -1)
    for (l_tmp = 0; l_tmp < nl ; l_tmp++) {
      scalar pp = ppl[l_tmp];
      pp[right]  = dirichlet(vpg[l_tmp]*x - upg[l_tmp]*y);
      pp[left]   = dirichlet(vpg[l_tmp]*x - upg[l_tmp]*y);
      pp[top]    = dirichlet(vpg[l_tmp]*x - upg[l_tmp]*y);
      pp[bottom] = dirichlet(vpg[l_tmp]*x - upg[l_tmp]*y);

      boundary({pp});
    }

}

/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

event init (i = 0)
{
  set_const();
}

/**
## Cleanup */

void trash_vars(){
  delete(pol), free (pol), pol = NULL;
  delete(qol), free (qol), qol = NULL;
#if MODE_PV_INVERT
  delete(pom), free (pom), pom = NULL;
  delete(qom), free (qom), qom = NULL;
  delete(cl2m), free (cl2m), cl2m = NULL;
  delete(cm2l), free (cm2l), cm2l = NULL;
  delete(iBul), free (iBul), iBul = NULL;
#endif
  delete(zetal), free (zetal), zetal = NULL;
  delete(qofl), free (qofl), qofl = NULL;
  delete(ppl), free (ppl), ppl = NULL;
  delete(Frl), free (Frl), Frl = NULL;
  delete(strl), free (strl), strl = NULL;
  delete(tmpl), free (tmpl), tmpl = NULL;
  free(dhf);
  free(dhc);
  free(idh0);
  free(idh1);
}

event cleanup (i = end, last) {
  trash_vars();
}

/**
   Python interface routines (should be in .i file but I need foreach)
 */

void pyset_field (scalar * psil, double * val1){
  int i = 0;
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar psi = psil[l];
      // transposed index
      int j = N*N*(i%nl) + N*((int)floor(i/nl)%N) + (int)floor(i/(N*nl));
      psi[] =  val1[j];
      i++;
    }
  boundary(psil);
}

void pyget_field (scalar * psil, double * val1){
  int i = 0;
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar psi = psil[l];
      // in principle, if the loop were reversed, it should be this line
//      int j = (i*nl*N)%(nl*N*N) + (nl*(int)floor(i/N))%(nl*N) + (int)floor(i/(N*N));
      // transposed index
      int j = N*N*l + N*((int)floor(i/nl)%N) + (int)floor(i/(N*nl));
      val1[j] = psi[];
      i++;
    }
}
