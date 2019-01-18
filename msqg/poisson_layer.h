/**
# 3d elliptic solver in the layer framework

We wish to solve the elliptic equation
$$
L(a) = \nabla\cdot (\nabla a) + \Gamma a = b\, ,
$$
with $\Gamma$ the vertical stretching operator
$$
\Gamma = \frac{\partial }{\partial z}\left(  S \frac{\partial a }{\partial z} \right)
$$

and with vertical BC
$$
\frac{\partial a}{\partial z} = 0
$$
at $z = 0$ and $z = -1$.

In the layer configuration, we assume a small number of vertical grid
points compared to number of horizontal grid points.

We store the stretching coefficients $S$ in a scalar field and we do
not use a multigrid solver in the vertical but a simple Jacobi
relaxation.

This function does not work with embeded BC yet.
This function does not work with TREES yet.
 */


struct Poisson_layer {
  scalar * a;
  scalar * b;
  (const) face vector alpha;
  (const) scalar lambda;
  (const) scalar * str0l;
  (const) scalar * str1l;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  bool (* embed_flux) (Point, scalar, vector, bool, double *);
#endif
};




static void relax_layer (scalar * al, scalar * bl, int l, void * data)
{
  struct Poisson_layer * p = (struct Poisson_layer *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  (const) scalar * str0l = p->str0l;
  (const) scalar * str1l = p->str1l;

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

  int nl = list_len (al);

// JACOBI and embeded not implemented yet (only direct reuse)
  
  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. */
  
  foreach_level_or_leaf (l) {

    double t0[nl], t1[nl], t2[nl], rhs[nl];


    if (nl > 1){
      // upper layer
      int ll = 0;

      scalar a = al[ll];
      scalar b = bl[ll];

      scalar str1 = str1l[ll];

      rhs[ll] = - sq(Delta)*b[];
      t1[ll] = sq(Delta)*str1[];
      t2[ll] = -sq(Delta)*str1[];
      foreach_dimension() {
        rhs[ll] += alpha.x[1]*a[1] + alpha.x[]*a[-1];
        t1[ll] += alpha.x[1] + alpha.x[];
      }
      
      // intermediate layers
      for (int ll = 1; ll < nl-1 ; ll++) {
        
      scalar a = al[ll];
      scalar b = bl[ll];

      scalar str0 = str0l[ll];
      scalar str1 = str1l[ll];

      rhs[ll] = - sq(Delta)*b[];
      t1[ll] = sq(Delta)*(str0[] + str1[]);
      t2[ll] = -sq(Delta)*str1[];
      t0[ll] = -sq(Delta)*str0[];

      foreach_dimension() {
        rhs[ll] += alpha.x[1]*a[1] + alpha.x[]*a[-1];
        t1[ll] += alpha.x[1] + alpha.x[];
      }
      
      }
      // lower layer
      ll = nl-1;

      a = al[ll];
      b = bl[ll];

      scalar str0 = str0l[ll];

      rhs[ll] = - sq(Delta)*b[];
      t1[ll] = sq(Delta)*str0[];
      t0[ll] = -sq(Delta)*str0[];
      foreach_dimension() {
        rhs[ll] += alpha.x[1]*a[1] + alpha.x[]*a[-1];
        t1[ll] += alpha.x[1] + alpha.x[];
      }

    /**
    We can now solve the tridiagonal system using the [Thomas
    algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (ll = 1; ll < nl; ll++) {
      t1[ll] -= t0[ll]*t2[ll-1]/t1[ll-1];
      rhs[ll] -= t0[ll]*rhs[ll-1]/t1[ll-1];
    }
    a = al[nl-1];
    a[] = t0[nl-1] = rhs[nl-1]/t1[nl-1];
    for (ll = nl - 2; ll >= 0; ll--) {
      a = al[ll];
      a[] = t0[ll] = (rhs[ll] - t2[ll]*t0[ll+1])/t1[ll];
    }

    }
  }
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double residual_layer (scalar * al, scalar * bl, scalar * resl, void * data)
{
  struct Poisson_layer * p = (struct Poisson_layer *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  (const) scalar * str0l = p->str0l;
  (const) scalar * str1l = p->str1l;
  double maxres = 0.;


  int nl = list_len (al);

// TREE not implemented yet

#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*face_gradient_x (a, 0);
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {

    // upper layer
    int l = 0;
    
    scalar a1 = al[l];
    scalar a2 = al[l+1];
    scalar b = bl[l];
    scalar res = resl[l];
    
    scalar str1 = str1l[l];
    
    
    res[] = b[] - str1[]*a2[] + str1[]*a1[];
    foreach_dimension()
      res[] += (alpha.x[0]*face_gradient_x (a1, 0) -
		alpha.x[1]*face_gradient_x (a1, 1))/Delta;  
    
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
    
    // intermediate layers
    for (int l = 1; l < nl-1 ; l++) {
      
      scalar a0 = al[l-1];
      scalar a1 = al[l];
      scalar a2 = al[l+1];
      scalar b = bl[l];
      scalar res = resl[l];
      
      scalar str0 = str0l[l];
      scalar str1 = str1l[l];
      
      
      res[] = b[] - str0[]*a0[] - str1[]*a2[] + (str0[] + str1[])*a1[] ;
      foreach_dimension()
        res[] += (alpha.x[0]*face_gradient_x (a1, 0) -
                  alpha.x[1]*face_gradient_x (a1, 1))/Delta;  
      
      if (fabs (res[]) > maxres)
        maxres = fabs (res[]);
      
    }
    // lower layer
    l = nl-1;
    
    scalar a0 = al[l-1];
    a1 = al[l];
    b = bl[l];
    res = resl[l];
    
    scalar str0 = str0l[l];
    
    res[] = b[] - str0[]*a0[] + str0[]*a1[];
    foreach_dimension()
      res[] += (alpha.x[0]*face_gradient_x (a1, 0) -
		alpha.x[1]*face_gradient_x (a1, 1))/Delta;  
    
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);

#endif // !TREE    
#if EMBED
    if (p->embed_flux) {
      double c;
      if (p->embed_flux (point, a, alpha, false, &c))
	res[] += c;
      else
	a[] = c, res[] = 0.;
    }
#endif // EMBED    
    /* if (fabs (res[]) > maxres) */
    /*   maxres = fabs (res[]); */
  }
  boundary (resl);
  return maxres;
}




mgstats poisson_layer (struct Poisson_layer p)
{

  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction ({alpha,lambda});

  scalar * str0l = p.str0l;
  scalar * str1l = p.str1l;
  restriction (str0l);
  restriction (str1l);

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar * a = p.a;
  scalar * b = p.b;
  mgstats s = mg_solve (a, b, residual_layer, relax_layer,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
