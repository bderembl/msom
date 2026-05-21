/**
# Multigrid Poisson--Helmholtz solvers

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

## Multigrid cycle

Here we implement the multigrid cycle proper. Given an initial guess
*a*, a residual *res*, a correction field *da* and a relaxation
function *relax*, we will provide an improved guess at the end of the
cycle. */

#include "vertex-utils.h"

trace
void mg_cycle (scalar * a, scalar * res, scalar * da,
	       void (* relax) (scalar * da, scalar * res,
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel, int maxlevel)
{


  boundary(res);
  /**
  We first define the residual on all levels. */
  for (scalar s in res)
    s.restriction = restriction_coarsen_vert;

  restriction (res);
  for (int l = 0; l <= depth(); l++) {
    boundary_level(res, l);
  }

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_vertex_level (l)
	for (scalar s in da)
	  foreach_blockf (s)
	    s[] = 0.;

    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */
    else {
      boundary_level (da, l - 1);
      foreach_vertex_level (l)
        for (scalar s in da)
          foreach_blockf (s)
            s[] = bilinear_vertex (point, s);
    }
    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    for (int i = 0; i < nrelax; i++) {
      boundary_level (da, l);
      relax (da, res, l, data);
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */
  boundary_level (da, maxlevel);

  foreach_vertex() {
    vertex scalar s, ds;
    for (s, ds in a, da) 
      foreach_blockf (s)
	s[] += ds[];
  }
  boundary(a); //vertex not automatic
}

/**
## Multigrid solver

The multigrid solver itself uses successive calls to the multigrid
cycle to refine an initial guess until a specified tolerance is
reached.

The maximum number of iterations is controlled by *NITERMAX* and the
tolerance by *TOLERANCE* with the default values below. */

int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3 [*];

/**
Information about the convergence of the solver is returned in a structure. */

typedef struct {
  int i;              // number of iterations
  double resb, resa;  // maximum residual before and after the iterations
  double sum;         // sum of r.h.s.
  int nrelax;         // number of relaxations
  int minlevel;       // minimum level of the multigrid hierarchy
} mgstats;

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. The optional number of relaxations is *nrelax* and *res* is
an optional list of fields used to store the residuals. The minimum
level of the hierarchy can be set (default is zero i.e. the root
cell). */

trace
mgstats mg_solve (scalar * a, scalar * b,
		  double (* residual) (scalar * a, scalar * b, scalar * res,
				       void * data),
		  void (* relax) (scalar * da, scalar * res, int depth,
				  void * data),
		  void * data = NULL,
		  int nrelax = 4,
		  scalar * res = NULL,
		  int minlevel = 1, // on vertex minlevel = 1
		  double tolerance = TOLERANCE)
{
  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  vertex scalar * da = list_clone (a), * pres = res;
  if (!res)
    res = list_clone (b);


  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];

  /* for (vertex scalar s in da){ */
  /*   s.restriction = restriction_vert; */
  /*   s.prolongation = refine_vert; */
  /* } */

  /* for (vertex scalar s in res){ */
  /*   s.prolongation = refine_vert; */
  /*   // BD: In principle, there is no need to define a BC for the rhs (or */
  /*   // residual).  However, for the Multigrid method, these BC are necessary to */
  /*   // converge in case a and b do not have the same BC (e.g. for the QG model with no slip BC) */
  /*   s[left]   = 0.; */
  /*   s[right]  = 0.; */
  /*   s[top]    = 0.; */
  /*   s[bottom] = 0.; */
  /* } */


  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0};
  double sum = 0.;
  vertex scalar rhs = b[0];
  foreach_vertex (reduction(+:sum))
    sum += rhs[];
  s.sum = sum;
  s.nrelax = nrelax > 0 ? nrelax : 4;

  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = (* residual) (a, b, res, data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > tolerance);
       s.i++) {
    mg_cycle (a, res, da, relax, data,
	      s.nrelax,
	      minlevel,
	      grid->maxdepth);
    s.resa = (* residual) (a, b, res, data);
    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa > tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	s.nrelax--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb = s.resa;
  }
  s.minlevel = minlevel;

  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > tolerance) {
    scalar v = a[0]; // fixme: should not be necessary
    fprintf (ferr,
	     "src/poisson.h:%d: warning: convergence for %s not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d tolerance: %g\n", LINENO, v.name,
	     s.i, s.resa, s.sum, s.nrelax, tolerance), fflush (ferr);
  }

  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!pres)
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
initial number of relaxations (default is one), *minlevel* controls
the minimum level of the hierarchy (default is one) and *res* is an
optional list of fields used to store the final residual (which can be
useful to monitor convergence). */

struct Poisson {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  double (* embed_flux) (Point, scalar, vector, double *);
#endif
};


