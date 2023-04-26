/**

For the barotropic solver
$$
q = \nabla^2 \psi
$$


*/

#include "wavelet_vertex.h"


vertex scalar psi_pg;
vertex scalar S2;
vertex scalar zeta;
vertex scalar psi_f;

vertex scalar topo[];
vertex scalar sig_lev[];



double idh0[nl_max] = {1.};
double idh1[nl_max] = {1.};
double N2[nl_max] = {1.};

int nbar = 0;

double scale_topo = 1.;

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
#define L_filt (Lfmax + (y/L0)*(Lfmin - Lfmax))

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
    point.l = 0;
    stretch[] = add*stretch[] + fac*S2[]*(psi[0,0,1] - psi[])*idh1[point.l];
    
    // intermediate layers
    for (point.l = 1; point.l < nl-1 ; point.l++) {
      stretch[] = add*stretch[] + fac*(S2[0,0,-1]*(psi[0,0,-1] - psi[])*idh0[point.l] + S2[]*(psi[0,0,1] - psi[])*idh1[point.l]);
    }
    
    // lower layer
    point.l = nl-1;
    stretch[] = add*stretch[] + fac*S2[0,0,-1]*(psi[0,0,-1] - psi[])*idh0[point.l];
    point.l = 0;
  }
  
  boundary({stretch});
}


trace
void rhs_pv_baroclinic(scalar q, scalar psi, scalar dqdt)
{
  
  foreach_vertex() {
    foreach_layer(){
      q[] *= mask[];
      psi[] *= mask[];
    }
  }


  comp_del2(psi, zeta, 0., 1.);

  foreach_vertex(){

    // upper layer
    point.l = 0;

    double ju,jd;

    jd = jacobian_l1(psi, psi) + jacobian_l1(psi_pg, psi) + jacobian_l1(psi, psi_pg);

    dqdt[] = -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - S2[]*jd*idh1[point.l]  \
      - beta_effect(psi);

      // intermediate layers
      for (point.l = 1; point.l < nl-1 ; point.l++) {

        ju = -jd;
        jd = jacobian_l1(psi, psi) + jacobian_l1(psi_pg, psi) + jacobian_l1(psi, psi_pg);
        
        dqdt[] = -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - S2[]*jd*idh1[point.l] - S2[0,0,-1]*ju*idh0[point.l] \
          - beta_effect(psi);

      }

      // lower layer
      point.l = nl-1;
      ju = -jd;
      dqdt[] = -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - S2[0,0,-1]*ju*idh0[point.l]  \
        - beta_effect(psi);

      /**
         Bottom friction and topography
      */
        dqdt[] += - hEkb*f_var/(2*dh[nl-1])*zeta[] - jacobian(psi, topo)*f_var/dh[nl-1];


      point.l = 0;
  }

  /**
     Dissipation
   */
  comp_stretch(zeta, dqdt, 1., nu);
  comp_del2(zeta, dqdt, 1., nu);

  /**
     Surface forcing
   */
  foreach_vertex() {
    dqdt[] += q_forcing[];
  }


  foreach_vertex() {
    foreach_layer(){
      dqdt[] *= mask[];
    }
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
    point.l = 0;

    rhs[point.l] = - sq(Delta)*b[]*mask[];
    t2[point.l] = -sq(Delta)*S2[]*idh1[point.l]*mask[];
    t1[point.l] = -t2[point.l];
    foreach_dimension() {
      rhs[point.l] += (a[1] + a[-1])*mask[];
      t1[point.l] += 2;
    }

    // intermediate layers
    for (point.l = 1; point.l < nl-1 ; point.l++) {
      rhs[point.l] = - sq(Delta)*b[]*mask[];
      t0[point.l] = -sq(Delta)*S2[0,0,-1]*idh0[point.l]*mask[];
      t2[point.l] = -sq(Delta)*S2[]*idh1[point.l]*mask[];
      t1[point.l] = -t0[point.l] - t2[point.l];

      foreach_dimension() {
        rhs[point.l] += (a[1] + a[-1])*mask[];
        t1[point.l] += 2;
      }
    }

    // lower layer
    point.l = nl-1;

    rhs[point.l] = - sq(Delta)*b[]*mask[];
    t0[point.l] = -sq(Delta)*S2[0,0,-1]*idh0[point.l];
    t1[point.l] = -t0[point.l];
    foreach_dimension() {
      rhs[point.l] += (a[1] + a[-1])*mask[];
      t1[point.l] += 2;
    }

    /**
       We can now solve the tridiagonal system using the [Thomas
       algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (point.l = 1; point.l < nl; point.l++) {
      t1[point.l] -= t0[point.l]*t2[point.l-1]/t1[point.l-1];
      rhs[point.l] -= t0[point.l]*rhs[point.l-1]/t1[point.l-1];
    }

    // lower layer
    point.l = nl-1;
    a[] = t0[nl-1] = rhs[nl-1]/t1[nl-1];
    for (point.l = nl - 2; point.l >= 0; point.l--) {
      a[] = t0[point.l] = (rhs[point.l] - t2[point.l]*t0[point.l+1])/t1[point.l];
    }
  point.l = 0;
  }
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
    point.l = 0;

      res[] = (b[] + S2[]*(a[] - a[0,0,1])*idh1[point.l])*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
      }
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);

    // intermediate layers
    for (point.l = 1; point.l < nl-1 ; point.l++) {
      res[] = (b[] + S2[0,0,-1]*(a[] - a[0,0,-1])*idh0[point.l] - S2[]*(a[0,0,1] - a[])*idh1[point.l])*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
      }
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);
    }

    // lower layer
    point.l = nl-1;

    res[] = (b[] + S2[0,0,-1]*(a[] - a[0,0,-1])*idh0[point.l])*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
      }
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);
      point.l = 0;
    }

    return maxres;

}

/**
   TODO: refactor wavelet filter
 */

trace
void wavelet_filter(scalar q, scalar psi)
{

  invert_q(psi,q);

  foreach_layer(){
    scalar w[];
    scalar psi_i[];

    w[top]    = 0;
    w[bottom] = 0;
    w[right]  = 0;
    w[left]   = 0;

    psi_i[top]    = dirichlet(0.);
    psi_i[bottom] = dirichlet(0.);
    psi_i[right]  = dirichlet(0.);
    psi_i[left]   = dirichlet(0.);



    // TODO: temporary fix to do wavelet transform with centered field
    // should be done in vertex field
    foreach()
      psi_i[] = 0.25*(psi[] + psi[1] + psi[0,1] + psi[1,1]);
    
    wavelet_mask(psi_i,w);
//    wavelet(psi,w);

    for (int l = 0; l <= depth(); l++) {
      foreach_vertex_level (l)
        w[] *= sig_lev[];
      boundary_level ({w}, l);
    }
    inverse_wavelet_mask (psi_i, w);
//    inverse_wavelet (psi, w);

    if (Lfmax < HUGE)
      foreach_vertex(){
        double psi_loc = 0.25*(psi_i[] + psi_i[-1] + psi_i[0,-1] + psi_i[-1,-1]);
        psi_f[] = (psi_f[]*nbar + psi_loc/dtflt)/(nbar+1);
        psi[] = (psi[] - psi_loc)*mask[];
      }

    /* foreach_vertex() */
    /*   psi[] *= mask[]; */
    
  }
  
  boundary({psi});

  comp_q(psi,q);
  nbar++;

}

/**
   Filter
*/
event filter (t = dtflt; t <= tend+1e-10;  t += dtflt) {
  fprintf(stdout,"Filter solution\n");
  wavelet_filter ( q, psi);
}

event defaults (i = 0){

  S2 = new vertex scalar[nl];
  psi_pg = new vertex scalar[nl];
  zeta = new vertex scalar[nl];
  psi_f = new vertex scalar[nl];

  S2.restriction = restriction_vert;
  S2.prolongation = refine_vert;

  psi_pg.restriction = restriction_vert;
  psi_pg.prolongation = refine_vert;

  zeta.restriction = restriction_vert;
  zeta.prolongation = refine_vert;

  psi_pg[left]   = psi_bc;
  psi_pg[right]  = psi_bc;
  psi_pg[bottom] = psi_bc;
  psi_pg[top]    = psi_bc;

  S2[left]   = neumann(0.);
  S2[right]  = neumann(0.);
  S2[bottom] = neumann(0.);
  S2[top]    = neumann(0.);


  reset ({S2, psi_pg, zeta, psi_f}, 0.);
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

  // defined in params.in
  foreach_vertex(){
    for (point.l = 0; point.l < nl-1 ; point.l++)
      S2[] = N2[point.l];
    point.l = 0;
    }

  fprintf(stdout, "Read input files:\n");

  FILE * fp;
  char name[80];
  sprintf (name,"input_vars_%dl_N%d.nc", nl,N);
  if ((fp = fopen(name, "r"))) {
    char N2_name[80] = "N2"; // trick to read N2 instead of S2
    char * N2_sav = S2.name; 
    S2.name = N2_name;
    read_nc({S2, psi_pg, mask, topo}, name);
//    read_nc({S2, psi_pg}, name);
    S2.name = N2_sav;
    fclose(fp);
    backup_file(name);
    fprintf(stdout, "%s .. ok\n", name);
  }


  /**
     Replace value S2: N^2 -> f^2/N^2
   */

  foreach_vertex(){
    for (point.l = 0; point.l < nl-1 ; point.l++)
      S2[] = sq(f_var)/S2[];
    point.l = 0;
  }

  restriction({S2});
  for (int l = 0; l <= depth(); l++) {
    boundary_level({S2}, l);
  }


  /**
     Rescale topo (TODO: temporary)
   */

  foreach_vertex()
    topo[] *= scale_topo;

  /**
     Compute filter

   */

  // low pass filter
  for (int l = depth(); l >= 0; l--) {
    foreach_vertex_level (l) {
      double ref_flag = 0;
      if (l < depth())
        foreach_child()
          ref_flag += sig_lev[];
      if (ref_flag > 0)
        sig_lev[] = 1;
      else{

        double L_filt2 = 0;
        if (fac_filt_Rd > 0){
          L_filt2 =  fac_filt_Rd*dh[0]/sqrt(S2[]);
        }
        else
          L_filt2 = L_filt;

        if (L_filt2 > 2*Delta)
          sig_lev[] = 0;
        else if (L_filt2 <= 2*Delta && L_filt2 > Delta)
          sig_lev[] = 1-(L_filt2-Delta)/Delta;
        else
          sig_lev[] = 1;
      }
    }
    boundary_level ({sig_lev}, l);
  }

  /* // high pass filter */
  /* for (int l = depth(); l >= 0; l--) { */
  /*   foreach_vertex_level (l) */
  /*     sig_lev[] = 1 - sig_lev[]; */
  /*   boundary_level ({sig_lev}, l);   */
  /* } */



    mask_c[top]    = dirichlet(0.);
    mask_c[bottom] = dirichlet(0.);
    mask_c[right]  = dirichlet(0.);
    mask_c[left]   = dirichlet(0.);



    // TODO: temporary fix to do wavelet transform with centered field
    // should be done in vertex field
    foreach()
      mask_c[] = 0.25*(mask[] + mask[1] + mask[0,1] + mask[1,1]);

    restriction({mask_c});
    for (int l = 0; l <= depth(); l++) {
      boundary_level({mask_c}, l);
    }

}


event cleanup (t = end, last)
{
  delete ({psi_pg, S2, zeta, psi_f});
}


event defaults (i = 0){

  rhs_pv = rhs_pv_baroclinic;
  comp_q = comp_q_baroclinic;
  invert_q = invert_q_baroclinic;
  relax_nodal = relax_baroclinic;
  residual_nodal = residual_baroclinic;
}
