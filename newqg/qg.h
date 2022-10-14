

/**
Quasi-Geostrophic code for Basilisk

This version is for a single layer but will eventually become a multi layer model

Solves


$$
\partial_t q + \nabla (u q ) + \beta v =  - \frac{h_{EK}*f0}{h_n} q  + nu \nabla^2 q + F
$$
with $q$ the vorticity
$$
q = \nabla^2 \psi
$$
$\psi$ the stream function
$$
u = -\frac{\partial \psi}{\partial y}
$$
and 
$$
v = \frac{\partial \psi}{\partial x}
$$

$h_{EK}$ is the thickness of the bottom Ekman layer, $H$ the thickness of the
lowermost layer and $f0$ is the Coriolis parameter. $\nu$ is the harmonic
viscosity.

The forcing $F$ is for now hard coded to
$$
F = \frac{3 \tau_0 \pi}{2 h_1 L} \sin(2*\pi*y/L0)*sin(\pi*y/L0);
$$
which correspond to a wind stress curl that will force a double gyre
configuration (cyclonic at the north, anticyclonic at the south).

The advective term $\nabla (uq)$ is coded with the Arakawa Jacobian which
conserves energy and enstrophy (see Cushman Roisin and Beckers 2011)


All parameters can be changed via a "params.in" file which should be in the
location as the executable. exemple of params.in:

N  = 64
L0 = 100

tau0 = 1e-3
nu = 5
beta = 0.5
dh   = [1.0]
sbc = 0

# timestepping
DT    = 5.e-2
tend  = 200.
dtout = 10
CFL   = 0.2
TOLERANCE = 1e-5


sbc is a parametre that controls the type of boundary conditions:
sbc = -1  -> periodic BC (experimental)
sbc = 0   -> free slip BC
sbc = big number (>10) -> no slip BC

For this single layer version: $q$ and $zeta$ are the same field (not true for
multiple layers)

*/

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
double gp_low = 0.;

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


// even if iRd2_low is a 2d field, it has to be defined as a layer field for the
// poisson solver
scalar iRd2_low;

/**
   Define useful operators:

   laplacian(psi) = \nabla^2 psi

   jacobian(psi, q) = dpsi/dx.dq/dy - dpsi/dy.dq/dx

   beta_effect(psi) = beta dpsi/dx = beta v
*/

#define laplacian(p) (p[1] + p[-1] + p[0,1] + p[0,-1] - 4*p[])/(sq(Delta))

// Arakawa Jacobian
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

// beta v
#define beta_effect(p) (beta*(p[1] - p[-1])/(2*Delta))

/** 
    Dynamical core
*/

trace
void invertq(scalar psi, scalar q)
{
/**
   Invert the poisson equation $\nabla^2 \psi - \psi/Rd^2 = q$
*/

  mgpsi = poisson(psi, q, lambda=iRd2_low);
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
void comp_stretch(scalar psi, scalar stretch, double add, double fac)
{

  if (gp_low != 0.) {
    // WARNING: iRd2_low only defined in first layer
    foreach()
      stretch[0,0,nl-1] = add*stretch[0,0,nl-1] + fac*(iRd2_low[]*psi[0,0,nl-1]);
  }

  boundary({stretch});
}

trace
void comp_q(scalar psi, scalar q)
{
  comp_del2  (psi, q, 0., 1.);
  comp_stretch(psi, q, 1., 1.);
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

/**
  Viscous dissipation
*/

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
   Prototype function for PV forcing term
*/
void surface_forcing  (scalar dqdt);


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

  iRd2_low = new scalar[nl]; 

  if (sbc >= 0) {
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
  }

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
  // Warning iRd2_low defined on upper layer only
  // TODO: to be upddated for multi layer
  if (gp_low != 0.) {
    foreach()
      iRd2_low[] = -sq(f0)/(gp_low*dh[nl-1]);
  } else {
    foreach()
      iRd2_low[] = 0.;
  }

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
  delete ({iRd2_low});
  free(evolving);
//  free(dhc);

}

event cleanup (t = end, last)
{
  trash_vars();
}
