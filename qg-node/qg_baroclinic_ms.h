/**

For the barotropic solver
$$
q = \nabla^2 \psi
$$


*/

vertex scalar psi_pg;
vertex scalar S2;
//vertex scalar N2; // alias for S2 (used to load variable only)
vertex scalar zeta;

vertex scalar topo[];


double idh0[nl_max] = {1.};
double idh1[nl_max] = {1.};


/**
   Jacobian (p[l], q[l+1])
*/

#define jacobian_l1(p,q) ((( p[1, 0   ]-p[-1, 0])*(q[0, 1,1]-q[ 0 ,-1,1]) \
                           +(p[0 ,-1  ]-p[ 0 ,1])*(q[1, 0,1]-q[-1, 0 ,1]) \
                           + p[1, 0   ]*( q[1,1 ,1] - q[1,-1 ,1])         \
                           - p[-1, 0  ]*( q[-1,1,1] - q[-1,-1,1])         \
                           - p[ 0 ,1  ]*( q[1,1 ,1] - q[-1,1 ,1])         \
                           + p[0 ,-1  ]*( q[1,-1,1] - q[-1,-1,1])         \
                           + q[ 0 ,1,1]*( p[1,1 ] - p[-1,1 ])         \
                           - q[0 ,-1,1]*( p[1,-1] - p[-1,-1])         \
                           - q[1, 0 ,1]*( p[1,1 ] - p[1,-1 ])         \
                           + q[-1, 0,1]*( p[-1,1] - p[-1,-1]))        \
                          /(12.*Delta*Delta))


#define f_var (f0 + flag_ms*beta*(y-0.5*L0))

/**
   Boundary condition for relative vorticity (same as PV)
   set_bc_ms will call the default set_bc
*/

void set_bc_ms()
{

  set_bc();

  zeta[left]   = 2*bc_fac/sq(Delta)*(psi[1]   - psi_bc);
  zeta[right]  = 2*bc_fac/sq(Delta)*(psi[]    - psi_bc);
  zeta[bottom] = 2*bc_fac/sq(Delta)*(psi[0,1] - psi_bc);
  zeta[top]    = 2*bc_fac/sq(Delta)*(psi[]    - psi_bc);

}


/**
   Stretching term
 */

trace
void comp_stretch(scalar psi, scalar stretch, double add, double fac)
{

  foreach_vertex() {

    // upper layer
    _layer = 0;
    stretch[] = add*stretch[] + fac*S2[]*(psi[0,0,1] - psi[])*idh1[_layer];
    
    // intermediate layers
    for (_layer = 1; _layer < nl-1 ; _layer++) {
      stretch[] = add*stretch[] + fac*(S2[0,0,-1]*(psi[0,0,-1] - psi[])*idh0[_layer] + S2[]*(psi[0,0,1] - psi[])*idh1[_layer]);
    }
    
    // lower layer
    _layer = nl-1;
    stretch[] = add*stretch[] + fac*S2[0,0,-1]*(psi[0,0,-1] - psi[])*idh0[_layer];
  } _layer = 0;
  
  boundary({stretch});
}


trace
void rhs_pv_baroclinic(scalar q, scalar psi, scalar dqdt)
{
  
  comp_del2(psi, zeta, 0., 1.);

  foreach_vertex(){

    // upper layer
    _layer = 0;

    double ju,jd;

    jd = jacobian_l1(psi, psi) + jacobian_l1(psi_pg, psi) + jacobian_l1(psi, psi_pg);

    dqdt[] = -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - S2[]*jd*idh1[_layer]  \
      - beta_effect(psi);

      // intermediate layers
      for (_layer = 1; _layer < nl-1 ; _layer++) {

        ju = -jd;
        jd = jacobian_l1(psi, psi) + jacobian_l1(psi_pg, psi) + jacobian_l1(psi, psi_pg);
        
        dqdt[] = -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - S2[]*jd*idh1[_layer] - S2[0,0,-1]*ju*idh0[_layer] \
          - beta_effect(psi);
      }

      // lower layer
      _layer = nl-1;
      ju = -jd;
      dqdt[] = -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - S2[0,0,-1]*ju*idh0[_layer]  \
        - beta_effect(psi);
      
  } _layer = 0;


  /**
     Dissipation
   */
  comp_stretch(zeta, dqdt, 1., nu);
  comp_del2(zeta, dqdt, 1., nu);

  /**
     Bottom friction and topography
   */
  foreach_vertex() {
    dqdt[0,0,nl-1] += - hEkb*f0/(2*dh[nl-1])*zeta[0,0,nl-1]\
      - jacobian(psi, topo)*f_var/dh[nl-1];
  }

  /**
     Surface forcing
   */
  foreach_vertex() {
    dqdt[] += q_forcing[];
  }



//    dqdt[] += -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - beta_effect(psi);

}


trace
void comp_q_baroclinic(scalar psi, scalar q)
{
  foreach_vertex()
    foreach_layer()
      q[] = laplacian(psi);

  comp_stretch(psi, q, 1., 1.);

  // need to set bc to force boundary conditions
  set_bc_ms();
  boundary({q});
}

/**
   Invert the poisson equation $\nabla^2 \psi = q$
*/

trace
void invert_q_baroclinic(scalar psi, scalar q)
{
  mgpsi = vpoisson(psi, q);
  // need to reset the values of the BC (has to do with vertices?)
  set_bc_ms();
  boundary({psi, q});
}


static void relax_baroclinic (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;

  foreach_vertex_level(l) {

    double t0[nl], t1[nl], t2[nl], rhs[nl];

    // upper layer
    _layer = 0;

    rhs[_layer] = - sq(Delta)*b[]*mask[];
    t2[_layer] = -sq(Delta)*S2[]*idh1[_layer]*mask[];
    t1[_layer] = -t2[_layer];
    foreach_dimension() {
      rhs[_layer] += (a[1] + a[-1])*mask[];
      t1[_layer] += 2;
    }

    // intermediate layers
    for (_layer = 1; _layer < nl-1 ; _layer++) {
      rhs[_layer] = - sq(Delta)*b[]*mask[];
      t0[_layer] = -sq(Delta)*S2[0,0,-1]*idh0[_layer]*mask[];
      t2[_layer] = -sq(Delta)*S2[]*idh1[_layer]*mask[];
      t1[_layer] = -t0[_layer] - t2[_layer];

      foreach_dimension() {
        rhs[_layer] += (a[1] + a[-1])*mask[];
        t1[_layer] += 2;
      }
    }

    // lower layer
    _layer = nl-1;

    rhs[_layer] = - sq(Delta)*b[]*mask[];
    t0[_layer] = -sq(Delta)*S2[0,0,-1]*idh0[_layer];
    t1[_layer] = -t0[_layer];
    foreach_dimension() {
      rhs[_layer] += (a[1] + a[-1])*mask[];
      t1[_layer] += 2;
    }

    /**
       We can now solve the tridiagonal system using the [Thomas
       algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (_layer = 1; _layer < nl; _layer++) {
      t1[_layer] -= t0[_layer]*t2[_layer-1]/t1[_layer-1];
      rhs[_layer] -= t0[_layer]*rhs[_layer-1]/t1[_layer-1];
    }

    // lower layer
    _layer = nl-1;
    a[] = t0[nl-1] = rhs[nl-1]/t1[nl-1];
    for (_layer = nl - 2; _layer >= 0; _layer--) {
      a[] = t0[_layer] = (rhs[_layer] - t2[_layer]*t0[_layer+1])/t1[_layer];
    }
  }
  _layer = 0;
  boundary_level({a}, l);

}

static double residual_baroclinic (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;

    foreach_vertex(reduction (max:maxres)) {

    // upper layer
    _layer = 0;

      res[] = (b[] + S2[]*(a[] - a[0,0,1])*idh1[_layer])*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
      }
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);

    // intermediate layers
    for (_layer = 1; _layer < nl-1 ; _layer++) {
      res[] = (b[] + S2[0,0,-1]*(a[] - a[0,0,-1])*idh0[_layer] - S2[]*(a[0,0,1] - a[])*idh1[_layer])*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
      }
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);
    }

    // lower layer
    _layer = nl-1;

    res[] = (b[] + S2[0,0,-1]*(a[] - a[0,0,-1])*idh0[_layer])*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
      }
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);

    } _layer = 0;

    return maxres;

}

event defaults (i = 0){

  S2 = new vertex scalar[nl];
  psi_pg = new vertex scalar[nl];
  zeta = new vertex scalar[nl];

  S2.restriction = restriction_vert;
  S2.prolongation = refine_vert;

  psi_pg.restriction = restriction_vert;
  psi_pg.prolongation = refine_vert;

  zeta.restriction = restriction_vert;
  zeta.prolongation = refine_vert;

  reset ({S2, psi_pg, zeta}, 0.);
}


event init (i = 0){

  /**
     Compute layer metrics
   */

  double dhc[nl_max] = {1.};

  for (int l = 0; l < nl-1; l++)
    dhc[l] = 0.5*(dh[l] + dh[l+1]);

  idh0[0] = 0.;
  idh1[0] = 1./(dhc[0]*dh[0]);
  for (int l = 1; l < nl-1 ; l++) {
    idh0[l] = 1./(dhc[l-1]*dh[l]);
    idh1[l] = 1./(dhc[l]*dh[l]);
  }
  idh0[nl-1] = 1./(dhc[nl-2]*dh[nl-1]);
  idh1[nl-1] = 0.;

  /**
     Read values of N^2
   */
  fprintf(stdout, "Read input files:\n");

  FILE * fp;
  char name[80];
  sprintf (name,"pgvars_%dl_N%d.nc", nl,N);
  if ((fp = fopen(name, "r"))) {
//    S2.name = "N2";
    read_nc({S2, psi_pg}, name);
//    S2.name = "S2";
    fclose(fp);
//    backup_file(name);
    fprintf(stdout, "%s .. ok\n", name);
  }

  /**
     Replace value S2: N^2 -> f^2/N^2
   */

  foreach_vertex()
    for (int l = 0; l < nl-1 ; l++) {
      S2[0,0,l] = sq(f_var)/S2[0,0,l];
    }
  restriction({S2});

}


event cleanup (t = end, last)
{
  delete ({psi_pg, S2});
}


event defaults (i = 0){

  flag_ms = 1;

  rhs_pv = rhs_pv_baroclinic;
  comp_q = comp_q_baroclinic;
  invert_q = invert_q_baroclinic;
  relax_nodal = relax_baroclinic;
  residual_nodal = residual_baroclinic;
}
