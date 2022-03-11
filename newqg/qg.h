#define LAYERS 1
#define nl_max 1000

double * dhf;
//double * dhc;

#include "predictor-corrector.h"
#include "poisson.h"

/**
   User defined constants
 */

double f0 = 1.0;
double beta = 0.;
double hEkb = 0.;
double tau0 = 0.;
double nu = 0.;
double sbc = 0.;
double tend = 1.;
double dtout = 1.; 
double dh[nl_max] = {1.};

char dpath[80]; // name of output dir
/** 
    Fields
*/

scalar q, psi, zeta;

/**
   Internal variables
 */
double bc_fac = 0;
scalar * evolving = NULL;
mgstats mgpsi;


/**
   Define useful operators:

   laplacian(psi) = \nabla^2 psi

   jacobian(psi, q) = dpsi/dx.dq/dy - dpsi/dy.dq/dx

   beta_effect(psi) = beta dpsi/dx = beta v
*/

#define laplacian(p) (p[1] + p[-1] + p[0,1] + p[0,-1] - 4*p[])/(sq(Delta))

#define jacobian(p,q) ((( p[1, 0 ]-p[-1, 0])*(q[0, 1]-q[ 0 ,-1]) \
                        +(p[0 ,-1]-p[ 0 ,1])*(q[1, 0]-q[-1, 0 ]) \
                        + p[1, 0 ]*( q[1,1 ] - q[1,-1 ])         \
                        - p[-1, 0]*( q[-1,1] - q[-1,-1])         \
                        - p[ 0 ,1]*( q[1,1 ] - q[-1,1 ])         \
                        + p[0 ,-1]*( q[1,-1] - q[-1,-1])         \
                        + q[ 0 ,1]*( p[1,1 ] - p[-1,1 ])         \
                        - q[0 ,-1]*( p[1,-1] - p[-1,-1])         \
                        - q[1, 0 ]*( p[1,1 ] - p[1,-1 ])         \
                        + q[-1, 0]*( p[-1,1] - p[-1,-1]))        \
                       /(12.*Delta*Delta))

#define beta_effect(p) (beta*(p[1] - p[-1])/(2*Delta))

/** 
    Dynamical core
*/

trace
void invertq(scalar psi, scalar q)
{
  mgpsi = poisson(psi, q);
  boundary({psi});
  boundary({q});
}

trace
void comp_del2(scalar psi, scalar zeta, double add, double fac)
{
  foreach()
    foreach_layer()
      zeta[] = add*zeta[] + fac*laplacian(psi);
  
  boundary({zeta});
}

trace
void comp_q(scalar psi, scalar q)
{
  comp_del2  (psi, q, 0., 1.);
  /* comp_stretch(psi, q, 1., 1.); */
  boundary({q});
}


/**
   Advection
*/

trace
double advection_pv(scalar zeta, scalar q, scalar psi, scalar dqdt, double dtmax)
{
  foreach()
    dqdt[] += -jacobian(psi, zeta) - beta_effect(psi);

  // compute dtmax (ajusted from timestep.h)
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax)){
    foreach_layer() {
      double u = 0.25*(psi[0,1] - psi[0,-1] + psi[-1,1] - psi[-1,-1])/Delta;
      if (u != 0.) {
        double dt = Delta/fabs(u);
        if (dt < dtmax) dtmax = dt;
      }
    }
  }

  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}

trace
void dissip  (scalar zeta, scalar dqdt)
{
  comp_del2(zeta, dqdt, 1., nu);
}

/**
   Bottom  Ekman Friction
*/

trace
void ekman_friction  (scalar zeta, scalar dqdt)
{
  foreach()
    dqdt[] -= hEkb*f0/(2*dh[nl-1])*zeta[];
}

/**
   surface forcing
*/
trace
void surface_forcing  (scalar dqdt)
{
  foreach()
    dqdt[] -= tau0/dh[0]*3/2*pi/L0*sin(2*pi*y/L0)*sin(pi*y/L0);
}


static void advance_qg (scalar * output, scalar * input,
                        scalar * updates, double dt)
{
  scalar qi = input[0];
  scalar qo = output[0];
  scalar dq = updates[0];

  foreach()
    foreach_layer() {
    qo[] = qi[] + dq[]*dt;
  }
//  boundary(output);
}


double update_qg (scalar * evolving, scalar * updates, double dtmax)
{
  // TODO check order loop
  foreach()
    for (scalar s in updates)
      foreach_layer()
        s[] = 0.;

  scalar q = evolving[0];
  scalar dqdt = updates[0];

  invertq(psi, q);
  comp_del2(psi, zeta, 0., 1.0);
  dtmax = advection_pv(zeta, q, psi, dqdt, dtmax);

  dissip(zeta, dqdt);
  ekman_friction(zeta, dqdt);
  surface_forcing(dqdt);
  
  return dtmax;
}


/**
   Declaration, initialize
*/

void set_vars()
{

  // derived const
  bc_fac = sbc/((0.5*sbc + 1)*sq(L0/N)); //TODO: no slip and adaptative mesh are not compatible

  q = new scalar[nl];
  psi = new scalar[nl];
  zeta = new scalar[nl];

  psi[right]  = dirichlet(0);
  psi[left]   = dirichlet(0);
  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);

  /* // TODO: replace by "ghost" */
  zeta[right]  = bc_fac*(psi[] - psi[1]);
  zeta[left]   = bc_fac*(psi[] - psi[-1]);
  zeta[top]    = bc_fac*(psi[] - psi[0,1]);
  zeta[bottom] = bc_fac*(psi[] - psi[0,-1]);

  q[right]  = bc_fac*(psi[] - psi[1]);
  q[left]   = bc_fac*(psi[] - psi[-1]);
  q[top]    = bc_fac*(psi[] - psi[0,1]);
  q[bottom] = bc_fac*(psi[] - psi[0,-1]);


  reset ({q, psi, zeta}, 0.);

//  if (nl > 1)
//    dhc = malloc ((nl-1)*sizeof(double));


  /**
     We overload the default 'advance' and 'update' functions of the
     predictor-corrector scheme and (TODO) setup the prolongation and
     restriction methods on trees. */

  evolving = list_copy({q});
  advance = advance_qg;
  update = update_qg;
  fprintf(stdout,"ok\n");
}

event defaults (i = 0){
  set_vars();
}
/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

void set_const() {
  comp_q(psi,q); // last part of init: invert PV
  boundary (all);
}

event init (i = 0) {  
  set_const(); 
}

/**
## Cleanup */

void trash_vars(){
  delete ({q, psi, zeta});
  free(evolving);
//  free(dhc);

}

event cleanup (t = end, last)
{
  trash_vars();
}
