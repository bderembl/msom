/**
# Multigrid Poisson--Helmholtz solvers

BD Comment: This file is very silimar to poisson.h. we use a different
version in order to experiment different boundary conditions.

We want to solve Poisson--Helmholtz equations of the general form
$$
L(a) = \nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
This can be done efficiently using a multigrid solver. 

An important aspect of Poisson--Helmholtz equations is that the
operator $L()$ is linear. This property can be used to build better
estimates of a solution by successive *corrections* to an initial
guess. If we define an approximate solution $\tilde{a}$ as
$$
\tilde{a} + da = a
$$
where $a$ is the exact (unknown) solution, using the linearity of the
operator we find that $da$ verifies
$$
L(da) = b - L(\tilde{a})
$$
where the right-hand-side is often called the *residual* of the
approximate solution $\tilde{a}$.



/**


## Multigrid cycle

Here we implement the multigrid cycle proper. Given an initial guess
*a*, a residual *res*, a correction field *da* and a relaxation
function *relax*, we will provide an improved guess at the end of the
cycle. */

// user defined boundary conditions
void circ_bc(scalar * psil);

void mg_cycle (scalar * a, scalar * res, scalar * da,
	       void (* relax) (scalar * da, scalar * res, 
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel, int maxlevel)
{

  /**
  We first define the residual on all levels. */

  restriction (res);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_level_or_leaf (l)
	for (scalar s in da)
	  s[] = 0.;

    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */

    else
      foreach_level (l)
	for (scalar s in da)
	  s[] = bilinear (point, s);

    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      s[] += ds[];
  }
  circ_bc (a);
//  boundary (a);
}

/**
## Multigrid solver

The multigrid solver itself uses successive calls to the multigrid
cycle to refine an initial guess until a specified tolerance is
reached. 

The maximum number of iterations is controlled by *NITERMAX* and the
tolerance by *TOLERANCE* with the default values below. */

int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;
int MINLEVEL = 0;

/**
Information about the convergence of the solver is returned in a structure. */

typedef struct {
  int i;              // number of iterations
  double resb, resa;  // maximum residual before and after the iterations
  double sum;         // sum of r.h.s.
  int nrelax;         // number of relaxations
} mgstats;

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. The optional number of relaxations is *nrelax* (default is
one) and *res* is an optional list of fields used to store the
residuals. */

struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
		       void * data);
  void (* relax) (scalar * da, scalar * res, int depth, 
		  void * data);
  void * data;
  
  int nrelax;
  scalar * res;
};

mgstats mg_solve (struct MGSolve p)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    for (scalar s in p.a) {
      scalar r = new scalar;
      res = list_append (res, r);
    }

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > TOLERANCE);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data, s.nrelax,
	      MINLEVEL, grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */
    
    if (s.resa > TOLERANCE) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	s.nrelax--;
    }
    resb = s.resa;
  }

  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > TOLERANCE)
    fprintf (ferr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d\n", 
	     s.i, s.resa, s.sum, s.nrelax), fflush (ferr);

  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

/**
## Application to the Poisson--Helmholtz equation

We now apply the generic multigrid solver to the Poisson--Helmholtz equation
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
We first setup the data structure required to pass the extra
parameters $\alpha$ and $\lambda$. We define $\alpha$ as a face
vector field because we need values at the face locations
corresponding to the face gradients of field $a$. 

*alpha* and *lambda* are declared as *(const)* to indicate that the
function works also when *alpha* and *lambda* are constant vector
(resp. scalar) fields. If *tolerance* is set, it supersedes the
default *TOLERANCE* of the multigrid solver, *nrelax* controls the
initial number of relaxations (default is one) and *res* is an optional
list of fields used to store the final residual (which can be useful
to monitor convergence). */

struct Poisson {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax;
  scalar * res;
  double minlevel;
};

/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;

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
    double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
    foreach_dimension() {
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
      d += alpha.x[1] + alpha.x[];
    }
    c[] = n/d;
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

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1])/Delta;
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] += (g.x[] - g.x[1])/Delta;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] += ((alpha.x[1] + alpha.x[])*a[]
		- alpha.x[1]*a[1] - alpha.x[]*a[-1])/sq(Delta);
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif
  boundary (resl);
  return maxres;
}

/**
## User interface

Finally we provide a generic user interface for a Poisson--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$ */

mgstats poisson (struct Poisson p)
{

  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x.i) {
    const vector alpha[] = {1.,1.,1.};
    p.alpha = alpha;
  }
  if (!p.lambda.i) {
    const scalar lambda[] = 0.;
    p.lambda = lambda;
  }

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction ({alpha,lambda});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  if (p.minlevel)
    MINLEVEL = p.minlevel;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ({a}, {b}, residual, relax, &p, p.nrelax, p.res);

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}

/**
## Projection of a velocity field

The function below "projects" the velocity field *u* onto the space of
divergence-free velocity fields i.e.
$$
\mathbf{u}^{n+1} \leftarrow \mathbf{u} - \Delta t\alpha\nabla p
$$
so that
$$
\nabla\cdot\mathbf{u}^{n+1} = 0
$$
This gives the Poisson equation for the pressure
$$
\nabla\cdot(\alpha\nabla p) = \frac{\nabla\cdot\mathbf{u}_*}{\Delta t}
$$ */

struct Project {
  face vector u;
  scalar p;
  face vector alpha; // optional: default unityf
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};

trace
mgstats project (struct Project q)
{
  face vector u = q.u;
  scalar p = q.p;
  (const) face vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_*$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[1] - u.x[];
    div[] /= dt*Delta;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha,
			 tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_{n+1}$ using $\mathbf{u}_*$ and $p$. */

  foreach_face()
    u.x[] -= dt*alpha.x[]*(p[] - p[-1])/Delta;
  boundary ((scalar *){u});

  return mgp;
}
