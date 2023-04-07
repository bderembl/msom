/**
Quasi-Geostrophic code for Basilisk

This version is for a single layer but will eventually become a multi layer model

Solves


$$
\partial_t q + \nabla (u q ) + \beta v =  - \frac{h_{EK}*f0}{h_n} q  + nu \nabla^2 q + F
$$
with $q$ the vorticity, $\psi$ the stream function
$$
u = -\frac{\partial \psi}{\partial y}
$$
and 
$$
v = \frac{\partial \psi}{\partial x}
$$


The definition of the potential vorticity $q$ is 
$$
q = \nabla^2 \psi + \frac{\partial}{\partial z} \frac{f^2}{N^2} \frac{\partial q}{\partial z}
$$

but reduces to 

$$
q = \nabla^2 \psi
$$

in the barotropic case

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
bc_fac = 0

# timestepping
DT    = 5.e-2
tend  = 200.
dtout = 10
CFL   = 0.2
TOLERANCE = 1e-5


bc_fac is a parameter that controls the type of boundary conditions:
bc_fac = 0      -> free slip BC (default)
bc_fac = [0..1] -> partial slip
bc_fac = 1      -> no slip BC
bc_fac = -1     -> periodic BC (experimental)

dtout is the frequency for writing outputs
tend the total time of the simulation

This version uses a prototype poisson solver for nodal fields written by Antoon
Van Hoft.  Hence all fields are defined at cell vertices. (still experimental)

*/


/**
   Anticipate layer version and so define generic foreach_all loop that can
   handle both the single layer and multi layer case
 */

//#define LAYERS 0
#define nl_max 1000


#include "predictor-corrector.h"
#include "nodal-poisson.h"


/**
   User defined constants
 */

double f0 = 1.0;
double beta = 0.;
double hEkb = 0.;
double tau0 = 0.;
double nu = 0.;
//double sbc = 0.; //retired
double tend = 100.;
double dtout = 1.; 
double dh[nl_max] = {1.};
double gp_l0 = 0.;

int flag_ms = 0;

char dpath[80]; // name of output dir
/** 
    Fields
*/

#if LAYERS
vertex scalar psi;
vertex scalar q;
#else
vertex scalar psi[];
vertex scalar q[];
#endif

vertex scalar q_forcing[]; 



/**
   Internal variables
 */
double bc_fac = 0;
double psi_bc = 0; // might change if we consider mass conservation
scalar * evolving = NULL;
mgstats mgpsi;


/**
   Energy diagnostics
*/

double dtdiag = -1; // non zero

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


void set_bc()
{
  // derived const
  //bc_fac = sbc/((0.5*sbc + 1)); // retired

  psi[left]   = psi_bc;
  psi[right]  = psi_bc;
  psi[bottom] = psi_bc;
  psi[top]    = psi_bc;

  /**
     First interior point minus boundary point in vertex convention.
   */

  q[left]   = 2*bc_fac/sq(Delta)*(psi[1]   - psi_bc);
  q[right]  = 2*bc_fac/sq(Delta)*(psi[]    - psi_bc);
  q[bottom] = 2*bc_fac/sq(Delta)*(psi[0,1] - psi_bc);
  q[top]    = 2*bc_fac/sq(Delta)*(psi[]    - psi_bc);
}


trace
void comp_del2(scalar psi, scalar zeta, double add, double fac)
{
  foreach_vertex()
#if LAYERS
    foreach_layer()
#endif
      zeta[] = add*zeta[] + fac*laplacian(psi);
  
  boundary({zeta});
}
/**
   Prototype functions
 */

void (* invert_q) (scalar psi, scalar q) = NULL;

void (* comp_q) (scalar psi, scalar q) = NULL;

void (* rhs_pv) (scalar q, scalar psi, scalar dqdt) = NULL;


/**
  compute dtmax (ajusted from timestep.h)

*/
double adjust_dt(scalar psi, double dtmax){
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax)){
#if LAYERS
    foreach_layer() {
#endif
      double u = (psi[0,1] - psi[])/Delta;
      if (u != 0.) {
        double dt = Delta/fabs(u);
        if (dt < dtmax) dtmax = dt;
      }
#if LAYERS
    }
#endif
  }

  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
/**
  Advance and update functions for time stepping. 
  Advance is different from the generic function because the loop is over vertices.

*/

static void advance_qg (scalar * output, scalar * input,
                        scalar * updates, double dt)
{

  vertex scalar qi = input[0];
  vertex scalar qo = output[0];
  vertex scalar dq = updates[0];

  foreach_vertex()
#if LAYERS
    foreach_layer()
#endif
      qo[] = qi[] + dq[]*dt;
//  boundary(output);
}


double update_qg (scalar * evolving, scalar * updates, double dtmax)
{

  vertex scalar q = evolving[0];
  vertex scalar dqdt = updates[0];

  invert_q(psi, q);
  rhs_pv(q, psi, dqdt);  
  dtmax = adjust_dt(psi, dtmax);

  return dtmax;
}


/**
   Energy diagnostics
*/

event write_1d_diag (t=0; t <= tend+1e-10; t += dtdiag){
  if (dtdiag > 0){
    FILE * fp;
    char name[80];
    sprintf (name,"%sdiag_1d.dat", dpath);

    if (i == 0){
      // create file
      if (pid() == 0) {
        if ((fp = fopen(name, "a"))) {
          fprintf(fp, "# time, ke, dissipation, forcing\n");
          fclose(fp);
        }
      }
    }
    
/**
   Using a foreach loop because foreach_vertex will count subdomain boundary
   points twice in MPI. Make sure the field is zero at the boundary
 */
    double ke = 0;
    double d_ke = 0;
    double f_ke = 0;

    // TODO: update in layer
    foreach(reduction(+:ke) reduction(+:d_ke) reduction(+:f_ke)){
      ke -= 0.5*psi[]*laplacian(psi)*sq(Delta);
      d_ke -= nu*psi[]*laplacian(q)*sq(Delta);
      f_ke -= psi[]*q_forcing[]*sq(Delta);
    }

    if (pid() == 0) {
      if ((fp = fopen(name, "a"))) {
        fprintf(fp, "%e, %e, %e, %e\n", t, ke, d_ke, f_ke);
        fclose(fp);
      }
    }
  } // dtdiag > 0
}



void set_vars()
{

#if LAYERS
  q = new vertex scalar[nl];
  psi = new vertex scalar[nl];
#endif

  psi.restriction = restriction_vert;
  psi.prolongation = refine_vert;

  q.restriction = restriction_vert;
  q.prolongation = refine_vert;


  reset ({q, psi}, 0.);
  reset ({q_forcing}, 0.);

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
  set_bc();
  set_vars();
}

/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

void set_const() {


  // Warning iRd2_l defined on upper layer only
  // TODO: to be upddated for multi layer
  /* foreach_vertex() */
  /*   S2[] = 0; */
    
//      iRd2_l[] = gp_l[] != 0 ? -sq(f0 + flag_ms*beta*(y-0.5*L0))/(gp_l[]*dh[nl-1]) : 0;

  comp_q(psi,q); // last part of init: compute PV (initial condition in psi)

  boundary (all);
}

event init (i = 0) {
  set_const(); 
}

/**
## Cleanup */

void trash_vars(){
#if LAYERS
  delete ({q, psi});
#endif
  free(evolving);

}

event cleanup (t = end, last)
{
  trash_vars();
}
