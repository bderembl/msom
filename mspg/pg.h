#include "predictor-corrector.h"
#include "elliptic.h"
#include "timestep.h"

scalar * bl  = NULL;
vector * ul = NULL;
scalar * pl  = NULL;
scalar * b_mel  = NULL;
vector * u_mel = NULL;
scalar * b_forcl  = NULL;

scalar * evolving = NULL;

// vertical grid
int nl = 1;
double * sc;
double * sf;
double ds;

// physical parameters
double r = 0.1;
double a = 0.2;
double kd = 3e-4;
double ys = 0; // southern latitude

// forcing
scalar b_surf[];
double tau_s = 1e-2;
double tau0 = 0.12;

// timesteping
double tend = 1; // end time
double dtout = 1; // Delat T output 
int nme = 0;

// 2d variables
scalar gg[];
scalar psibt[];
scalar wind_effect[];
face vector ubt[];

// diffusion coefficient
face vector kf[];

// elliptic solver
mgstats mgD;
face vector ronh[];
vector fonh[];
double omega = 0.3;


char outdir[80] = "./";

// user defined diffusivity coeff
double k (double x, double y, double s);

// user defined wind and wind derivative
double taux   (double x, double y);
double tauy   (double x, double y);
double taux_y (double x, double y);
double tauy_x (double x, double y);


/**
## Circular boundary condition

   implement a new type of BC. For now, assume square domain
  */


void circ_bc(scalar * psil)
{

  /**  temporary: assume square basin */
  int nbc = 4*N;

  for (scalar psi in psil) {

  double ad[nbc], bd[nbc], cd[nbc], rhs[nbc], sol[nbc];
    
  int ibc = 0;
  foreach_boundary (bottom){
    rhs[ibc] = psi[];
    ad[ibc] = -y/(4*r);
    bd[ibc] = 1.0;
    cd[ibc] = -ad[ibc];    
    ibc++;
  }
  foreach_boundary (right){
    rhs[ibc] = psi[];
    ad[ibc] = -y/(4*r);
    bd[ibc] = 1.0;
    cd[ibc] = -ad[ibc];    
    ibc++;
  }
  ibc += N-1;
  foreach_boundary (top){
    rhs[ibc] = psi[];
    ad[ibc] = -y/(4*r);
    bd[ibc] = 1.0;
    cd[ibc] = -ad[ibc];    
    ibc--;
  }
  ibc += 2*N;
  foreach_boundary (left){
    rhs[ibc] = psi[];
    ad[ibc] = -y/(4*r);
    bd[ibc] = 1.0;
    cd[ibc] = -ad[ibc];    
    ibc--;
  }
    /**
       We can now solve the tridiagonal system using the [Thomas
       algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (int ibc = 1; ibc < nbc; ibc++) {
      bd[ibc] -= ad[ibc]*cd[ibc-1]/bd[ibc-1];
      rhs[ibc] -= ad[ibc]*rhs[ibc-1]/bd[ibc-1];
    }
    

    sol[nbc-1] = ad[nbc-1] = rhs[nbc-1]/bd[nbc-1];
    for (int ibc = nbc - 2; ibc >= 0; ibc--) {
      sol[ibc] = ad[ibc] = (rhs[ibc] - cd[ibc]*ad[ibc+1])/bd[ibc];
    }

  ibc = 0;
  foreach_boundary (bottom){
//    if (ibc == 0) { psi[-1,-1] = 0.; } // corner
    psi[ghost] = 2*sol[ibc] - psi[];
    ibc++;
  }
  foreach_boundary (right){
    psi[ghost] = 2*sol[ibc] - psi[];
//    if (ibc == N) { psi[1,-1] = 0.5*(psi[1,0] + psi[0,-1]); } // corner
    if (ibc == N) { psi[1,-1] = (psi[1,0] + psi[0,-1] - psi[]); } // corner
    ibc++;
  }
  ibc += N-1;
  foreach_boundary (top){
    psi[ghost] = 2*sol[ibc] - psi[];
    /* if (ibc == 2*N) { psi[1,1] += 0.5*(psi[1,0] + psi[0,1]); } // corner */
    if (ibc == 2*N) { psi[1,1] = (psi[1,0] + psi[0,1] - psi[]); } // corner
    ibc--;
  }
  ibc += 2*N;
  foreach_boundary (left){
    psi[ghost] = 2*sol[ibc] - psi[];
    if (ibc == 3*N) { psi[-1,1] = (psi[-1,0] + psi[0,1] - psi[]); } // corner
    if (ibc == 4*N-1) { psi[-1,-1] = (psi[-1,0] + psi[0,-1] - psi[]); } // corner
    /* if (ibc == 3*N) { psi[-1,1] += 0.5*(psi[-1,0] + psi[0,1]); } // corner */
    /* if (ibc == 4*N-1) { psi[-1,-1] += 0.5*(psi[-1,0] + psi[0,-1]); } // corner */
    ibc--;
  }
}
}

void add_circ_bc(scalar * psil, scalar psibt) {

  for (scalar psi in psil) {

  int nbc = 4*N;

  int ibc = 0;
  foreach_boundary (bottom){
    psi[ghost] += psibt[ghost];
    ibc++;
  }
  foreach_boundary (right){
    psi[ghost] += psibt[ghost];
    if (ibc == N) { psi[1,-1] += 0.5*(psibt[1,0] + psibt[0,-1]); } // corner
    ibc++;
  }
  ibc += N-1;
  foreach_boundary (top){
    psi[ghost] += psibt[ghost];
    if (ibc == 2*N) { psi[1,1] += 0.5*(psibt[1,0] + psibt[0,1]); } // corner
    ibc--;
  }
  ibc += 2*N;
  foreach_boundary (left){
    psi[ghost] += psibt[ghost];
    if (ibc == 3*N) { psi[-1,1] += 0.5*(psibt[-1,0] + psibt[0,1]); } // corner
    if (ibc == 4*N-1) { psi[-1,-1] += 0.5*(psibt[-1,0] + psibt[0,-1]); } // corner
    ibc--;
  }
}
}

/**
## Elliptic Solver

   The routines for the elliptic solver are almost identical to the one
   in poisson.h. The main difference is the pseudo SOR method with the
   relaxation factor $\omega$. We added this because if the frictional
   term is weak the matrix wee want to invert is no longer diagonally
   dominant.
*/

struct Btsolver {
  scalar a, b;
  (const) face vector alpha;
  (const)  vector beta;
  (const) double omega;
  double tolerance;
  int nrelax;
  scalar * res;
  double minlevel;
};

static double residual_bt (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Btsolver * p = data;
  (const) face vector alpha = p->alpha;
  (const) vector beta = p->beta;
  double maxres = 0.;
  struct { double x, y; } f = {-1.,1.};
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1])/Delta;
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[];
    foreach_dimension() {
      res[] += (g.x[] - g.x[1])/Delta;
      res[] += f.x*beta.y[]*0.5*(a[1] - a[-1])/Delta;
    }
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[];
    foreach_dimension(){
      res[] += ((alpha.x[1] + alpha.x[])*a[]
  		- alpha.x[1]*a[1] - alpha.x[]*a[-1])/sq(Delta);
      res[] += f.x*beta.y[]*0.5*(a[1] - a[-1])/Delta;
    }

    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif
  boundary (resl);
  return maxres;
}


static void relax_bt (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Btsolver * p = data;
  (const) face vector alpha = p->alpha;
  (const) vector beta = p->beta;
  (const) double omega = p->omega;
  struct { double x, y; } f = {-1.,1.};

  /**
     We use either Jacobi (under)relaxation or we directly reuse values
     as soon as they are updated. For Jacobi, we need to allocate space
     for the new field *c*. Jacobi is useful mostly as it gives results
     which are independent of the order in which the cells are
     traversed. This is not the case for the simple traversal, which
     means for example that results will depend on whether a tree or
     a multigrid is used (because cells will be traversed in a different
     order). The same comment applies to OpenMP or MPI parallelism. In
     practice however Jacobi convergence tends to be slower than simple
     reuse. */
  
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  
  /**
     We use the face values of $\alpha$ to weight the gradients of the
     5-points Laplacian operator. We get the relaxation function. */

  foreach_level_or_leaf (l) {

    double n = - sq(Delta)*b[], d = 0.;
    foreach_dimension() {
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
      d += alpha.x[1] + alpha.x[];
      n -= f.x*beta.y[]*0.5*(a[1] - a[-1])*Delta;
    }
    c[] = (1-omega)*c[] + omega*n/d;
  }

  /**
     For weighted Jacobi we under-relax by using a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
  
#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = a[];
  trash ({a});
  foreach_level_or_leaf (l)
    a[] = a1[];
#endif
}

trace
mgstats btsolver (struct Btsolver p)
{

  face vector alpha = p.alpha;
  vector beta = p.beta;
  restriction ((scalar *){alpha,beta});

  /**
     If *tolerance* is set it supersedes the default of the multigrid
     solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  if (p.minlevel)
    MINLEVEL = p.minlevel;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ({a}, {b}, residual_bt, relax_bt, &p, p.nrelax, p.res);

  /**
     We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}


// compute BT velocity
double bt_velocity(face vector ubt, scalar psibt)
{

  struct { double x, y; } f = {-1.,1.};
  foreach_face(){
    
    ubt.x[] = (-r*((psibt[] - psibt[-1])/Delta )
               + f.x*y*(0.25*(psibt[0,1] - psibt[0,-1] + psibt[-1,1] - psibt[-1,-1])/Delta
         		))/(sq(r) + sq(y));
  }
  boundary ((scalar *){ubt});
}

/**
##  3d advection    

   Vertical boundary conditions (ghost points) needed for advection
   and horizontal diffusion because of the metric term in sigma
   coordinates
*/

void vertbc  (scalar * bl)
{
  scalar b0 = bl[0];
  scalar b1 = bl[1];
  foreach()
    b0[] = b1[];
  boundary({b0});

  b0 = bl[nl+1];
  b1 = bl[nl];
  foreach()
    b0[] =  b1[];
  boundary({b0});
}



/**
   the advection is a 3d flux divergence: 
*/

// TODO: add upper Ekman layer

trace

double advection (scalar * bl, vector * ul, scalar * dbl, double dtmax)
{
  
  scalar w0[];
  scalar w1[];
  // lower BC
  foreach()
    w0[] = 0.;

  vertbc(bl);
  foreach() {
    for (int l = nl; l > 0 ; l--) {
      
      scalar b = bl[l];
      scalar b1 = bl[l-1];
      scalar b0 = bl[l+1];
      scalar db = dbl[l];
      face vector u = ul[l];
      
      w1[] = w0[] - (u.x[1] - u.x[] + u.y[0,1] - u.y[])*ds/Delta;

      db[] += ((b[] + b[-1,0])*u.x[] -
      	       (b[] + b[1,0] )*u.x[1,0] +
      	       (b[] + b[0,-1])*u.y[] -
      	       (b[] + b[0,1] )*u.y[0,1]
      	       )/(2.*Delta) +
      	      ((b[] + b0[])*w0[] -
      	       (b[] + b1[])*w1[])/(2.*ds);

      w0[] = w1[];
    }
  }

  for (int l = nl; l > 0 ; l--) {
    face vector u = ul[l];
    dtmax = timestep (u, dtmax);
  }

  return dtmax;

}

/**
## diffusion

   The horizontal diffusion is computed explicitely and the vertical
   diffusion is computed implicitly. The treatment of vertical boundaries
   (air-sea interface and sea floor) is handled with ghots points which
   are adjusted in vertbc 

*/
trace
void vdiff_implicit(scalar * bl, double dt)
{
  double ad[nl], bd[nl], cd[nl], rhs[nl];
  
  /**
     surface BC */
  scalar b0 = bl[1];
  foreach()
    b0[] += dt*2*kd*k(x,y,sf[0])/(sq(ds))*b_surf[];
  
  foreach() {
    for (int l = 0; l < nl; l++) {
      b0 = bl[l+1];
      rhs[l] = b0[];
    }
    
    double K0 = kd*k(x,y,sf[0]);
    double K1 = kd*k(x,y,sf[1]);
    
    // upper layer
    ad[0] = 0.;
    cd[0] = - dt*K1/sq(ds);
    bd[0] = 1 + dt*K1/sq(ds) + 2*dt*K0/sq(ds);
    
    for (int l = 1; l < nl-1; l++) {
      
      K0 = K1;
      K1 = kd*k(x,y,sf[l+1]);
      
      ad[l] = - dt*K0/sq(ds);
      cd[l] = - dt*K1/sq(ds);
      bd[l] = 1. - ad[l] - cd[l];
      
    }
    /* Lower layer */
    K0 = K1;
    K1 = 0.;
    ad[nl-1] = - dt*K0/sq(ds);
    cd[nl-1] = - dt*K1/sq(ds);
    bd[nl-1] = 1. - ad[nl-1] - cd[nl-1];
    
    /**
       We can now solve the tridiagonal system using the [Thomas
       algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (int l = 1; l < nl; l++) {
      bd[l] -= ad[l]*cd[l-1]/bd[l-1];
      rhs[l] -= ad[l]*rhs[l-1]/bd[l-1];
    }
    
    b0 = bl[nl];
    b0[] = ad[nl-1] = rhs[nl-1]/bd[nl-1];
    for (int l = nl - 2; l >= 0; l--) {
      b0 = bl[l+1];
      b0[] = ad[l] = (rhs[l] - cd[l]*ad[l+1])/bd[l];
    }
  }
}

/**
   Horizontal diffusion
*/

trace
void hdiffusion  (scalar * bl, scalar * dbl) 
{
  vertbc(bl);

  face vector hdiff[];
  
  for (int l = 1; l < nl+1 ; l++) {
    foreach_face() {
      scalar b = bl[l];
      hdiff.x[] = sq(a)*kd*k(x,y,sc[l])*((b[] - b[-1])/Delta);
    }
  /**
     no flux side boundary conditions (ok with BC on b)*/
  /* hdiff.n[right]  = 0.0; */
  /* hdiff.n[left]   = 0.0; */
  /* hdiff.n[top]    = 0.0; */
  /* hdiff.n[bottom] = 0.0; */
//    boundary ((scalar *){hdiff});


  /**
     and divergence of the flux*/
  foreach(){
      scalar db = dbl[l];
      db[] += (hdiff.x[1,0]- hdiff.x[] + hdiff.y[0,1] - hdiff.y[])/Delta;
    }
  }
}

/**
## Convection
Restore unstable stratification to neutral profile */

void convection(scalar * bl)
{
  foreach() {
    /**
       instability criterion: convective instability occurs
       if b[l] < b[l+1]. **ASSUME LAYERS EQUALLY SPACED** */
    int l = 1;
    scalar b0 = bl[l];
    scalar b1 = bl[l+1];

    while(  l < nl) {
      if( b1[] > b0[]) {
        /**
           layer equaly spaced: standard averaging */
        b1[] = 0.5*(b0[] + b1[]);
        b0[] = b1[];
      }
      l++;
      b0 = bl[l];
      b1 = bl[l+1];
    }
  }
}


/**
## surface forcing
 */

void forcing_implicit(scalar * bl, double dt)
{

  scalar b = bl[1];
  foreach() 
    b[] = (b_surf[]*dt + b[]*tau_s)/(dt + tau_s);
}

/**
## QG forcing
 */

void qg_forcing(scalar * dbl, scalar * b_forcl)
{
  foreach() 
    for (int l = 1; l < nl+1 ; l++) {
      scalar db = dbl[l];
      scalar bf = b_forcl[l];
      db[] += bf[];
    }
}

/**
## Momentum
 */

void momentum(scalar * bl, vector * ul, vector * dul)
{

  foreach() {
    scalar p = pl[1];
    scalar b = bl[1];
    p[] = - b[]*0.5*ds;
  }

  foreach()
    for (int l = 2; l < nl+1 ; l++) {
    scalar b = bl[l];
    scalar b0 = bl[l-1];
    scalar p = pl[l];
    scalar p0 = pl[l-1];
    p[] = p0[] - 0.5*(b0[] + b[])*ds;
  }
  boundary(pl);


  struct { double x, y; } f = {-1.,1.};
  foreach_face() {
    for (int l = 1; l < nl+1 ; l++) {
      scalar b = bl[l];
      scalar p = pl[l];
      face vector u = ul[l];
      face vector du = dul[l];

      du.x[] = -(p[] - p[-1])/Delta
        - f.x*y*0.25*(u.y[] + u.y[-1] + u.y[0,1] + u.y[-1,1])
        - r*u.x[];
    }
  }
}


void adjust_bt_velocity(vector * ul, double btfac)
{

  face vector u_me[];
  foreach_face(){
    u_me.x[] = 0.;
    for (int l = 1; l < nl+1 ; l++) {
      face vector u = ul[l];
      u_me.x[] += u.x[]*ds;
    }
  }

  foreach_face(){
    for (int l = 1; l < nl+1 ; l++) {
      face vector u = ul[l];
      u.x[] += btfac*ubt.x[] - u_me.x[];
    }
  }

  for (int l = 1; l < nl+1 ; l++) {
    face vector u = ul[l];
    boundary ((scalar *){u});
  }
}

/**
## time stepping routines 

   We use the predictor corrector implementation */


static void advance_pg (scalar * output, scalar * input, 
			scalar * updates, double dt)
{
  scalar * bul = (scalar *) &updates[0];
  vector * uul = (vector *) &updates[nl+2];

  scalar * bil = (scalar *) &input[0];
  vector * uil = (vector *) &input[nl+2];

  scalar * bol = (scalar *) &output[0];
  vector * uol = (vector *) &output[nl+2];

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {
      scalar bi = bil[l];
      scalar bo = bol[l];
      scalar db = bul[l];
      bo[] = bi[] + db[]*dt;
    }
  }

  forcing_implicit(bol,dt);
  vdiff_implicit(bol,dt);
  convection(bol);

  foreach_face() {
    for (int l = 1; l < nl+1 ; l++) {
      face vector ui = uil[l];
      face vector uo = uol[l];
      face vector du = uul[l];
      uo.x[] = ui.x[] + du.x[]*dt;
    }
  }

  adjust_bt_velocity(uol,1);

  boundary (bol);

  for (int l = 1; l < nl+1 ; l++) {
    face vector uo = uol[l];
    boundary ((scalar *){uo});
  }

}

double update_pg (scalar * evolving, scalar * updates, double dtmax)
{
  foreach()
    for (scalar s in updates)
      s[] = 0.;

  scalar * bl = (scalar *) &evolving[0];
  vector * ul = (vector *) &evolving[nl+2];

  scalar * dbl = (scalar *) &updates[0];
  vector * dul = (vector *) &updates[nl+2];

  dtmax = advection (bl, ul, dbl, dtmax);
  hdiffusion(bl, dbl);
  qg_forcing(dbl,b_forcl);
  momentum(bl, ul, dul);

  return dtmax;
}

void set_vars()
{
  assert (pl     == NULL);
  assert (bl     == NULL);
  assert (ul     == NULL);
  assert (nl > 0);


  // nl+ 2 layers
  for (int l = 0; l < nl+2; l++) {
    scalar b = new scalar;
    bl  = list_append (bl, b);

    scalar p = new scalar;
    pl  = list_append (pl, p);

    face vector u = new face vector;
    ul = vectors_append (ul, u);

    face vector um = new face vector;
    u_mel = vectors_append (u_mel, um);

  }

  evolving = list_concat(bl, (scalar *) ul);

  // init vertical grid
  ds = 1./nl;
  sc = malloc (nl*sizeof(double));
  sf = malloc ((nl+1)*sizeof(double));
  
  sf[0] = -1.0;
  for (int l = 1; l < nl+1; l++)
    sf[l] = sf[l-1] + ds;

  sc[0] = -1.0 + 0.5*ds;
  for (int l = 1; l < nl; l++)
    sc[l] = sc[l-1] + ds;


  // elliptic solver
  foreach(){
    fonh.x[] = 0.;
    fonh.y[] = -(sq(r)-sq(y))/sq(sq(r)+sq(y));
  }
  foreach_face(){
    ronh.x[] = r/(sq(r)+sq(y));
  }

  foreach(){
    wind_effect[] = tau0*taux_y(x,y);  /* Samelson */
  }

  b_mel = list_clone(bl);
  b_forcl = list_clone(bl);

  /**
     Initialize variables */
  foreach() {
    for (scalar b in bl)
      b[] = 0.0;
  }
  foreach_face() {
    for (vector u in ul)
      u.x[] = 0.0;
  }

  foreach() 
    for (scalar b in b_mel)
      b[] = 0.0;
  
  foreach() 
    for (scalar b in b_forcl)
      b[] = 0.0;
  
  foreach_face() {
    for (vector u in u_mel)
      u.x[] = 0.0;
  }

  // TODO: not sure why but I need these lines (even with circ_bc)
  psibt[right]  = dirichlet(0);
  psibt[left]   = dirichlet(0);
  psibt[top]    = dirichlet(0);
  psibt[bottom] = dirichlet(0);

  
  // TODO change u.t = 0
  for (vector u in ul) {
   u.n[right]  = 0.0;
   u.n[left]   = 0.0;
   u.n[top]    = 0.0;
   u.n[bottom] = 0.0;
  }

   ubt.n[right]  = 0.0;
   ubt.n[left]   = 0.0;
   ubt.n[top]    = 0.0;
   ubt.n[bottom] = 0.0;

  advance = advance_pg;
  update = update_pg;

}


event defaults (i = 0){
  set_vars();
}


event init (i = 0)
{
  mgD = btsolver(psibt,wind_effect,ronh,fonh,omega);
  bt_velocity(ubt, psibt);

  boundary (all);
}


/**
## Cleanup */
void trash_vars(){
  free (bl), bl = NULL;
  free (pl), pl = NULL;
  free (ul), ul = NULL;
  free (evolving);
  free (sc);
  free (sf);
  free (b_mel), b_mel = NULL;
  free (u_mel), u_mel = NULL;
}

event cleanup (i = end, last) {
  trash_vars();
}

////////////////: python specific routines /////////////////////////

scalar * pyupdates = NULL;


double dtconv = 1e-1; // only needed for the bifurcation solver; do not choose it too small


// continuation parameters, etc
int contpar = 0;
double forcing_mag = 1.0;

/**
    Explicit vertical diffusion routine for the bifuraction solver */

trace
void vdiff_explicit  (scalar * bl, scalar * dbl) {
  vertbc(bl);
  foreach(){
    for (int l = 1; l < nl+1 ; l++) {
      scalar b0 = bl[l];
      scalar b1 = bl[l+1];
      scalar b2 = bl[l-1];
      scalar db = dbl[l];

      db[] +=   (kd*k(x,y,sf[l-1])*(b2[] - b0[])
               - kd*k(x,y,sf[l])*(b0[] - b1[]))/sq(ds);
    }
  }
}

void convection_tend(scalar * bl, scalar * dbl, double dt)
{

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {
      scalar b = bl[l];
      scalar b_sav = b_mel[l];
      b_sav[] = 1.0*b[];
    }
  }

  convection(bl);

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {
      scalar b = bl[l];
      scalar b_sav = b_mel[l];
      scalar db = dbl[l];
      db[] += (b[]-b_sav[])/dt;
      b[] = b_sav[];
    }
  }
}

void forcing_explicit(scalar * bl, scalar *dbl)
{

  scalar b = bl[1];
  scalar db = dbl[1];
  foreach()
    db[] += (b_surf[] - b[])/tau_s;
}


/**
   Python interface routines (should be in .i file but I need foreach)
 */
/**
   Set and adjust the continuation parameter via the python interface
 */
void pyset_contpar(int pycontpar) { 
  contpar = pycontpar; 
}

void pyadjust_contpar(double contpar_val){
  if (contpar == 1) {
    forcing_mag = contpar_val;
    foreach() 
      b_surf[] = forcing_mag*6*cos(pi*(y-ys));
  }
}

void pyinit_const(int pynl){ 
  nl = pynl;
  a = sqrt(3.0e-2/k(0,0,0)); 
  r = 0.02; 
  tau_s = 3.0e-2;
  omega = 0.2;

  ys = 0.3;
  origin (0.0, ys);

}

void pyinit_last(){

  pyupdates = list_clone(evolving);

  double tau0 = 0.12;
  foreach()
    wind_effect[] = tau0*(2*pi*tau0*cos(2*(y-ys)*pi));

  foreach()
    b_surf[] = 6*cos(pi*(y-ys));

  mgD = btsolver(psibt,wind_effect,ronh,fonh,omega);
  bt_velocity(ubt, psibt);

  boundary (all);
}

void pyset_field (scalar * evolving, double * val1, int len1){
  
  scalar * bl = (scalar *) &evolving[0];
  vector * ul = (vector *) &evolving[nl+2];

  int i = 0;
  foreach()
    for (int l = 1; l < nl+1; l++) {
    scalar b = bl[l];
      b[] =  val1[i];
      i++;
    }
  boundary(bl);
  vertbc(bl);

  foreach_face(x)
    for (int l = 1; l < nl+1 ; l++) {
      face vector u = ul[l];
      u.x[] =  val1[i];
      i++;
    }

  foreach_face(y)
    for (int l = 1; l < nl+1 ; l++) {
      face vector u = ul[l];
      u.y[] =  val1[i];
      i++;
    }

  for (int l = 1; l < nl+1 ; l++) {
    face vector u = ul[l];
    boundary ((scalar *){u});
  }
}


void pyget_field (scalar * evolving, double * val1, int len1){
  
  scalar * bl = (scalar *) &evolving[0];
  vector * ul = (vector *) &evolving[nl+2];

  int i = 0;
  foreach()
    for (int l = 1; l < nl+1; l++) {
    scalar b = bl[l];
      val1[i] = b[];
      i++;
    }

  foreach_face(x)
    for (int l = 1; l < nl+1 ; l++) {
      face vector u = ul[l];
       val1[i] = u.x[];
      i++;
    }

  foreach_face(y)
    for (int l = 1; l < nl+1 ; l++) {
      face vector u = ul[l];
       val1[i] = u.y[];
      i++;
    }
}


void pystep ( double * val1, int len1,
              double * val2, int len2){
  
  foreach()
    for (scalar s in pyupdates)
      s[] = 0.;

  scalar * bl = (scalar *) &evolving[0];
  vector * ul = (vector *) &evolving[nl+2];

  scalar * dbl = (scalar *) &pyupdates[0];
  vector * dul = (vector *) &pyupdates[nl+2];

  pyset_field ( evolving, val1, len1);

  double dtloc = 1.0;

  adjust_bt_velocity(ul,1);

  dtloc = advection (bl, ul, dbl, dtloc);
  hdiffusion(bl, dbl);
//  qg_forcing(dbl,b_forcl);
  forcing_explicit(bl,dbl);
  vdiff_explicit(bl,dbl);
  convection_tend(bl,dbl,dtconv);

 momentum(bl, ul, dul);
 adjust_bt_velocity(dul,0);

  pyget_field ( pyupdates, val2, len2);
} 

void pytrash_vars (){
  
  trash_vars();
  free (pyupdates);
}
