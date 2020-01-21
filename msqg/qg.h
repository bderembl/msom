/**
   multiple scale quasi-geostrophic model 
*/

double * dhf;
double * dhc;
double * idh0;
double * idh1;

#include "layer.h"
#include "predictor-corrector.h"
#include "poisson.h"
#include "poisson_layer.h"
#include "timestep.h"

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
scalar Ro[];
scalar Rd[];
scalar sig_filt[];
scalar sig_lev[];
int lsmin, lsmax; 
double afilt = 10.;  // filter size = afilt*Rd
double Lfmax = 10.;  // max filter length scale

// temporary psi field with samem BC as p
scalar * tmpl = NULL;

scalar * evolving = NULL;

mgstats mgpsi;

int nl = 1;

double Re  = 0.; // reynolds number
double Re4 = 0.; // bihormonic reynolds number
double Re6 = 0.; // bihormonic reynolds number
double iRe = 0.0; // inverse reynolds number
double iRe4 = 0.0; // inverse bihormonic reynolds number
double iRe6 = 0.0; // inverse bihormonic reynolds number
double Ekb = 1e-2; // Ekman number (bottom)
double Eks = 1e-2; // Ekman number (surface)
double Rom = 1e-2; // Mean Rossby number (if negative: Ro = cte = -Rom)
double Frm = 1e-2; // Mean Froude number
double beta = 0.5;
double tau0 = 0.; // wind stress curl
double sbc = 0.; // slip BC: 0: free slip (OK), big: no slip (TO BE FINISHED)

double tend = 1; // end time
double dtflt = 1; // Delat T filtering
double dtout = 1; // Delat T output

int nbar = 0;
int ediag = -1;  // ediag = -1: no ediag, 0: psi*dqdt, 1: (psi+pg)*dqdt

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
        zeta[ghost] = add*zeta[ghost] 
          + fac*sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
      foreach_boundary(right) {
        zeta[ghost] = add*zeta[ghost] 
          + fac*sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
      foreach_boundary(top) {
        zeta[ghost] = add*zeta[ghost] 
          + fac*sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
      foreach_boundary(bottom) {
        zeta[ghost] = add*zeta[ghost] 
          + fac*sbc/((0.5*sbc + 1)*sq(Delta))*(po[]-po[ghost]);}
    }
  }

}

trace
void comp_stretch(scalar * pol, scalar * stretchl, double add)
{

  foreach() {

    if (nl > 1){
      // upper layer
      int l = 0;
      scalar po_1 = pol[l];
      scalar po_2 = pol[l+1];
      scalar stretch = stretchl[l];
      scalar s1 = strl[l];

      stretch[] = add*stretch[] + s1[]*(po_2[] - po_1[])*idh1[l];

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        scalar po_0 = pol[l-1];
        scalar po_1 = pol[l];
        scalar po_2 = pol[l+1];
        scalar stretch = stretchl[l];
        scalar s0 = strl[l-1];
        scalar s1 = strl[l];

        stretch[] = add*stretch[] + s0[]*(po_0[] - po_1[])*idh0[l] + s1[]*(po_2[] - po_1[])*idh1[l];
      }
      // lower layer
      l = nl-1;
      scalar po_0 = pol[l-1];
      po_1 = pol[l];
      stretch = stretchl[l];
      scalar s0 = strl[l-1];

      stretch[] = add*stretch[] + s0[]*(po_0[] - po_1[])*idh0[l];
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
    uf.x[] = f.x*(po[] - po[0,-1])/Delta;
  }
}

trace
double advection  (scalar * pol, scalar * dqol, double dtmax)
{

  comp_del2(pol, zetal, 0., 1.0);

  foreach() {
    double ju,jd;

    if (nl > 1){

      // upper layer
      int l = 0;
      scalar zo  = zetal[l];
      scalar po  = pol[l];
      scalar pp  = ppl[l];
      scalar dqo = dqol[l];
      scalar po2  = pol[l+1];
      scalar pp2  = ppl[l+1];
      scalar s1 = strl[l];

      jd = jacobian(po, po2) + jacobian(pp, po2) + jacobian(po, pp2);
      dqo[] += jacobian(po, zo) + jacobian(pp, zo) + beta_effect(po) + s1[]*jd*idh1[l];

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        zo  = zetal[l];
        po  = pol[l];
        pp  = ppl[l];
        dqo = dqol[l];
        po2  = pol[l+1];
        pp2  = ppl[l+1];
        scalar s0 = strl[l-1];
        scalar s1 = strl[l];

        ju = -jd;
        jd = jacobian(po, po2) + jacobian(pp, po2) + jacobian(po, pp2);

        dqo[] += jacobian(po, zo) + jacobian(pp, zo) + beta_effect(po) + s0[]*ju*idh0[l] + s1[]*jd*idh1[l];
      }

      // lower layer
      l = nl-1;

      zo  = zetal[l];
      po  = pol[l];
      pp  = ppl[l];
      dqo = dqol[l];
      scalar s0 = strl[l-1];

      ju = -jd;
      dqo[] += jacobian(po, zo) + jacobian(pp, zo) + beta_effect(po) + s0[]*ju*idh0[l];
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
  //    dtmax = DT;
   }
  return dtmax;
}

trace
void comp_q(scalar * pol, scalar * qol)
{
  comp_del2  (pol, qol, 0., 1.);
  comp_stretch(pol, qol, 1.);
  // TODO: not the rigght BC if partial or no slip
  boundary(qol);
}


trace
void dissip  (scalar * qol, scalar * dqol)
{
  comp_del2(qol, tmpl, 0., 1.);

  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar dqo = dqol[l];
      scalar q2 = tmpl[l];
      dqo[] += q2[]*iRe;
    }

  comp_del2(tmpl, zetal, 0., 1.);

  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar dqo = dqol[l];
      scalar q4 = zetal[l];
      dqo[] += q4[]*iRe4;
    }

  comp_del2(zetal, dqol, 1., iRe6);
}

/**
   Bottom and surface Ekman Friction
*/

trace
void ekman_friction  (scalar * pol, scalar * dqol)
{
  foreach() {
    scalar dqos = dqol[0];
    scalar pos = pol[0];

    scalar dqob = dqol[nl-1];
    scalar pob = pol[nl-1];
    dqos[] -= Eks*laplacian(pos);
    dqob[] -= Ekb*laplacian(pob);
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
    dqo[] -= tau0*sin(2*pi*y/L0);
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
  dtmax = advection(pol, updates, dtmax);
  dissip(evolving, updates);
  ekman_friction(pol, updates);
  surface_forcing(updates);

  return dtmax;
}

/**
   Filter
*/
event filter (t = dtflt; t <= tend+1e-10;  t += dtflt) {
  fprintf(stdout,"Filter solution\n");
  wavelet_filter ( qol, pol, qofl, dtflt, nbar)
}

/**********************************************************************
*                       End of dynamical core                         *
***********************************************************************/

/**
   Read input parameters
 */

void read_params()
{
  FILE * fp;
  if ((fp = fopen("params.in", "rt"))) {
    char tempbuff[100];
    char tmps1[16];
    char tmps2[16];

    while(fgets(tempbuff,100,fp)) {
      sscanf(tempbuff, "%15s = %15s", tmps1, tmps2);
      if      (strcmp(tmps1,"N")    ==0) { N     = atoi(tmps2); }
      else if (strcmp(tmps1,"nl")   ==0) { nl    = atoi(tmps2); }
      else if (strcmp(tmps1,"ediag")==0) { ediag = atoi(tmps2); }
      else if (strcmp(tmps1,"L0")   ==0) { L0    = atof(tmps2); }
      else if (strcmp(tmps1,"Rom")  ==0) { Rom   = atof(tmps2); }
      else if (strcmp(tmps1,"Ekb")  ==0) { Ekb   = atof(tmps2); }
      else if (strcmp(tmps1,"Eks")  ==0) { Eks   = atof(tmps2); }
      else if (strcmp(tmps1,"Re")   ==0) { Re    = atof(tmps2); }
      else if (strcmp(tmps1,"Re4")  ==0) { Re4   = atof(tmps2); }
      else if (strcmp(tmps1,"Re6")  ==0) { Re6   = atof(tmps2); }
      else if (strcmp(tmps1,"beta") ==0) { beta  = atof(tmps2); }
      else if (strcmp(tmps1,"afilt")==0) { afilt = atof(tmps2); }
      else if (strcmp(tmps1,"Lfmax")==0) { Lfmax = atof(tmps2); }
      else if (strcmp(tmps1,"DT")   ==0) { DT    = atof(tmps2); }
      else if (strcmp(tmps1,"tend") ==0) { tend  = atof(tmps2); }
      else if (strcmp(tmps1,"dtout")==0) { dtout = atof(tmps2); }
      else if (strcmp(tmps1,"dtflt")==0) { dtflt = atof(tmps2); }
      else if (strcmp(tmps1,"CFL")  ==0) { CFL   = atof(tmps2); }
    }
    fclose(fp);
  } else {
    fprintf(stdout, "file params.in not found\n");
    exit(0);
  }
  if (Re  == 0) iRe  = 0.; else iRe  =  1/Re;
  if (Re4 == 0) iRe4 = 0.; else iRe4 = -1/Re4;
  if (Re6 == 0) iRe6 = 0.; else iRe6 =  1/Re6;

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
  fprintf(stdout,"Create main variables .. ");
  pol   = create_layer_var(pol,nl,0);
  qol   = create_layer_var(qol,nl,0);
  ppl   = create_layer_var(ppl,nl,0);
  zetal = create_layer_var(zetal,nl,0);
  tmpl  = create_layer_var(tmpl,nl,0);
  Frl   = create_layer_var(Frl,nl,1);
  strl = create_layer_var(strl,nl,1);
  qofl  = create_layer_var(qofl,nl,0);
#if MODE_PV_INVERT
  pom = create_layer_var(pom,nl,0);
  qom = create_layer_var(qom,nl,0);
  iBul  = create_layer_var(iBul,nl,1);
#endif

  /**
     Mode to layer inversion matrices (dimesnion: $nl^2$)
  */
#if MODE_PV_INVERT
  int nl2 = nl*nl;
  cl2m = create_layer_var(cl2m,nl2,1);
  cm2l = create_layer_var(cm2l,nl2,1);
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
    dhf[l] = 1/nl;

  foreach()
    for (scalar Fr in Frl)
      Fr[] =  Frm; 

  foreach(){
    Ro[] = Rom; 
    Rd[] = 1.; 
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

double shape(double x,double sigma){ return (1-exp(-0.5*sq(x/sigma)));}


void set_const() {

/**
   Layer thickness and large scale variables
*/
  fprintf(stdout, "Read input files:\n");

  char name[80];
  sprintf (name,"dh_%dl.bin", nl);
  float dh[nl];
  FILE * fp = fopen (name, "r");
  fread(&dh, sizeof(float), nl, fp);
  fclose(fp);

  for (int l = 0; l < nl ; l++)
    dhf[l] = dh[l];

  fprintf(stdout, "%s .. ok\n", name);

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


  foreach()
    Ro[] = (Rom > 0) ? Rom/(1 + Rom*beta*(y-0.5*L0)) : -Rom;

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
  // sponge
  /* foreach() */
  /*   sig_filt[] = min(Lfmax*shape(x,5.)*shape(x-L0,5.)*shape(y,5.0)*shape(y-L0,5.0),Lfmax); */
  /* foreach() */
  /*   sig_filt[] = min(Lfmax*shape(x,10.)*shape(x-L0,10.)*shape(y,10.0)*shape(y-L0,10.0),sig_filt[]); */
    

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

  /**
     BC for all fields */
  boundary (all);
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
  free (pol), pol = NULL;
  free (qol), qol = NULL;
#if MODE_PV_INVERT
  free (pom), pom = NULL;
  free (qom), qom = NULL;
  free (cl2m), cl2m = NULL;
  free (cm2l), cm2l = NULL;
  free (iBul), iBul = NULL;
#endif
  free (zetal), zetal = NULL;
  free (qofl), qofl = NULL;
  free (ppl), ppl = NULL;
  free (Frl), Frl = NULL;
  free (strl), strl = NULL;
  free (tmpl), tmpl = NULL;
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
