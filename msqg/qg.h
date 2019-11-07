/**
   multiple scale quasi-geostrophic model 
*/

#include "predictor-corrector.h"
#include "poisson.h"
#include "poisson_layer.h"
#include "timestep.h"

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
scalar * str0l = NULL;
scalar * str1l = NULL;
// If PV inversion is done with the tri-diagonal method, 
// then no need to refine cl2m,cm2l,iBul, Frl, Ro
scalar * cl2m  = NULL;
scalar * cm2l  = NULL;
scalar * iBul = NULL; // inverse burger number
scalar * Frl = NULL;
scalar Ro[];
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
double * dhf;
double * dhc;

double Re = 1e1; // reynolds number
double Re4 = 1e3; // bihormonic reynolds number
double Ek = 1e-2; // Ekman number
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

#include "eigmode.h"

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
  mgpsi = poisson_layer (pol, qol, str0l = str0l, str1l = str1l, tolerance = 1e-3);
#endif

  boundary(pol);
}

trace
void comp_del2(scalar * pol, scalar * zetal, double add, double fac)
{
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar po = pol[l];
      scalar zeta = zetal[l];
      zeta[] = add*zeta[] 
        + fac*(po[1] + po[-1] + po[0,1] + po[0,-1] - 4*po[])/(sq(Delta));
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
      scalar s1 = str1l[l];

      stretch[] = add*stretch[] + s1[]*po_2[] - s1[]*po_1[];

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        scalar po_0 = pol[l-1];
        scalar po_1 = pol[l];
        scalar po_2 = pol[l+1];
        scalar stretch = stretchl[l];
        scalar s0 = str0l[l];
        scalar s1 = str1l[l];

        stretch[] = add*stretch[] + s0[]*po_0[] + s1[]*po_2[] - (s0[] + s1[])*po_1[];
      }
      // lower layer
      l = nl-1;
      scalar po_0 = pol[l-1];
      po_1 = pol[l];
      stretch = stretchl[l];
      scalar s0 = str0l[l];

      stretch[] = add*stretch[] + s0[]*po_0[] - s0[]*po_1[];
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
double advection  (scalar * qol, scalar * pol, scalar * dqol, double dtmax)
{

  foreach() {
    double ju,jd;

    if (nl > 1){

      // upper layer
      int l = 0;
      scalar qo  = qol[l];
      scalar po  = pol[l];
      scalar pp  = ppl[l];
      scalar dqo = dqol[l];
      scalar po2  = pol[l+1];
      scalar pp2  = ppl[l+1];
      scalar s1 = str1l[l];

      jd = jacobian(po, po2) + jacobian(pp, po2) + jacobian(po, pp2);
      dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s1[]*jd;

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        qo  = qol[l];
        po  = pol[l];
        pp  = ppl[l];
        dqo = dqol[l];
        po2  = pol[l+1];
        pp2  = ppl[l+1];
        scalar s0 = str0l[l];
        scalar s1 = str1l[l];

        ju = -jd;
        jd = jacobian(po, po2) + jacobian(pp, po2) + jacobian(po, pp2);

        dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju + s1[]*jd;
      }

      // lower layer
      l = nl-1;

      qo  = qol[l];
      po  = pol[l];
      pp  = ppl[l];
      dqo = dqol[l];
      scalar s0 = str0l[l];

      ju = -jd;
      dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju;
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
void dissip  (scalar * zetal, scalar * dqol)
{
  double iRe = 1/Re;
  comp_del2(zetal, tmpl, 0., 1.);

  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar dqo = dqol[l];
      scalar p4 = tmpl[l];
      dqo[] += p4[]*iRe;
    }

  double iRe4 = -1/Re4;
  comp_del2(tmpl, dqol, 1., iRe4);
}

/**
   Ekman Friction on bottom layer
*/

trace
void bottom_friction  (scalar * zetal, scalar * dqol)
{
  scalar dqo = dqol[nl-1];
  scalar zeta = zetal[nl-1];
  foreach()
    dqo[] -= Ek*zeta[];
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


void list_copy_deep (scalar * listin, scalar * listout)
{
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar listi = listin[l];
      scalar listo = listout[l];
      listo[] = listi[];
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
  if (dtflt == -1.0){
    list_copy_deep (tmpl, qol);
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
  dtmax = advection(zetal, pol, updates, dtmax);
  dissip(zetal, updates);
  bottom_friction(zetal, updates);
  surface_forcing(updates);

  return dtmax;
}

/**  
   ## Layerered variables initialization, etc
*/

scalar * create_layer_var (scalar * psil, int nl, int bc_type)
{
  assert (psil == NULL);
  assert (nl > 0);

  for (int l = 0; l < nl; l++) {
    scalar po = new scalar;
    psil = list_append (psil, po);
    
    if (bc_type == 0){
      po[right]  = dirichlet(0);
      po[left]   = dirichlet(0);
      po[top]    = dirichlet(0);
      po[bottom] = dirichlet(0);
    }
  }

  foreach() 
    for (scalar po in psil) {po[] = 0.0;} 

  return psil;
}

void reset_layer_var(scalar *psil)
{
  foreach()
    for (scalar po in psil) {po[] = 0.0;}
}

void set_vars()
{
  pol   = create_layer_var(pol,nl,0);
  qol   = create_layer_var(qol,nl,0);
  ppl   = create_layer_var(ppl,nl,0);
  zetal = create_layer_var(zetal,nl,0);
  tmpl  = create_layer_var(tmpl,nl,0);
  Frl   = create_layer_var(Frl,nl,1);
  iBul  = create_layer_var(iBul,nl,1);
  str0l = create_layer_var(str0l,nl,1);
  str1l = create_layer_var(str1l,nl,1);
  qofl  = create_layer_var(qofl,nl,0);
//  create_layer_var(qosl,nl);
#if MODE_PV_INVERT
  pom = create_layer_var(pom,nl,0);
  qom = create_layer_var(qom,nl,0);
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
  for (int l = 0; l < nl; l++)
    dhf[l] = 1/nl;

  foreach()
    for (scalar Fr in Frl)
      Fr[] =  Frm; 

  foreach()
    Ro[] = Rom; 

  /**
     We overload the default 'advance' and 'update' functions of the
     predictor-corrector scheme and (TODO) setup the prolongation and
     restriction methods on trees. */

  advance = advance_qg;
  update = update_qg;
}

event defaults (i = 0){
  fprintf(stdout,"Create main variables .. ");
  set_vars();
  fprintf(stdout,"ok\n");
}


/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

event init (i = 0)
{
  for (int l = 0; l < nl-1; l++)
    dhc[l] = 0.5*(dhf[l] + dhf[l+1]);

  /* foreach() */
  /*   Ro[] = uref/((fref + betad*(y-0.5*L0)*lref)*lref); */
  foreach()
    Ro[] = (Rom > 0) ? Rom/(1 + Rom*beta*(y-0.5*L0)) : -Rom;

  /* foreach()  */
  /*   for (int l = 0; l < nl ; l++) { */
  /*     scalar pp  = ppl[l]; */
  /*     pp[] = pp[]*Ro[]; */
  /*   } */

  /**
     Compute vertical stretching coef matrix
   */

  foreach() {

    if (nl > 1){
      // upper layer
      int l = 0;
      scalar Fr1 = Frl[l];
      scalar str1 = str1l[l];
      str1[] = sq(Fr1[]/Ro[])/( dhc[l]*dhf[l]);

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
        
        scalar Fr0 = Frl[l-1];
        scalar Fr1 = Frl[l];
        
        scalar str0 = str0l[l];
        scalar str1 = str1l[l];

        str0[] = sq(Fr0[]/Ro[])/( dhc[l-1]*dhf[l]);
        str1[] = sq(Fr1[]/Ro[])/( dhc[l]*dhf[l]);

      }
      // lower layer
      l = nl-1;
      scalar Fr0 = Frl[l-1];
      scalar str0 = str0l[l];
      str0[] = sq(Fr0[]/Ro[])/( dhc[l-1]*dhf[l]);
    }
  }

  fprintf(stdout,"Compute vertical modes .. ");
  eigmod(dhf, dhc, Ro, Frl, cl2m, cm2l, iBul);
  fprintf(stdout,"ok\n ");

  /**
     compute filter length scale and wavelet coeffs*/
  scalar iRd = iBul[1];
  foreach()
    sig_filt[] = min(afilt*sqrt(-1/iRd[]),Lfmax);


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

  /**
     BC for all fields */
  boundary (all);
}


void write_field(scalar *psil, char name[], double rescale){

  int nl1 = list_len (psil);

  if (rescale != 0)
    foreach()
      for (int l = 0; l < nl1 ; l++) {
        scalar psi  = psil[l];
        psi[] *= rescale;
      }
  
  FILE * fp = fopen (name, "w");
  output_matrixl (psil, fp);
  fclose(fp);
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
#endif
  free (zetal), zetal = NULL;
  free (qofl), qofl = NULL;
  free (ppl), ppl = NULL;
  free (iBul), iBul = NULL;
  free (Frl), Frl = NULL;
  free (str0l), str0l = NULL;
  free (str1l), str1l = NULL;
  free (tmpl), tmpl = NULL;
//  free (qosl), qosl = NULL;
  free(dhf);
  free(dhc);
}

event cleanup (i = end, last) {
  trash_vars();
}

/**
   Energy diagnostic routines
*/
#include "qg_energy.h"
