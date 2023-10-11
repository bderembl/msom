/**

For the barotropic solver
$$
q = \nabla^2 \psi
$$


*/

// temporary topo
//vertex scalar topo[];


trace
void rhs_pv_barotropic(scalar q, scalar psi, scalar dqdt)
{

  foreach_vertex(){
    dqdt[] = -jacobian(psi, q)                     \
      - beta_effect(psi)                           \
      - hEkb*f0/(2*dh[nl-1])*q[]                   \
      + q_forcing[]                                \
//      - jacobian(psi, topo)*f0/dh[nl-1]            \
      + nu*laplacian(q);
  }
//    dqdt[] += -jacobian(psi, zeta) - jacobian(psi_pg, zeta) - beta_effect(psi);

}


trace
void comp_q_barotropic(scalar psi, scalar q)
{
  foreach_vertex()
    q[] = laplacian(psi);

  boundary({q});
}

/**
   Invert the poisson equation $\nabla^2 \psi = q$
*/

trace
void invert_q_barotropic(scalar psi, scalar q)
{

  mgpsi = vpoisson(psi, q);
//  mgpsi = vpoisson(psi, q, lambda=iRd2_l);
  // need to reset the values of the BC (has to do with vertices?)
  set_bc();
  boundary({psi, q});
}


static void relax_barotropic (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;

        foreach_vertex_level(l) {
          double d = 0;
//          double d = - lambda[]*sq(Delta);
          a[] = -b[]*sq(Delta);
          foreach_dimension() {
              a[] += (a[1] + a[-1])*mask[];
              d += 2.;
          }
          a[] /= d;
        }
        boundary_level({a}, l);
}

static double residual_barotropic (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;

    foreach_vertex(reduction (max:maxres)) {
      res[] = b[]*mask[];
//      res[] = (b[] - lambda[]*a[])*mask[];
      foreach_dimension() {
          res[] -= (a[-1] - 2.*a[] + a[1])/(sq(Delta))*mask[];
      }
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);
    }
    return maxres;

}



event defaults (i = 0){

  rhs_pv = rhs_pv_barotropic;
  comp_q = comp_q_barotropic;
  invert_q = invert_q_barotropic;
  relax_nodal = relax_barotropic;
  residual_nodal = residual_barotropic;

}



// temporary topo
/* event init (i = 0) { */

/*   FILE * fp; */
/*   char name[80]; */
/*   sprintf (name,"input_vars_%dl_N%d.nc", nl,N); */
/*   if ((fp = fopen(name, "r"))) { */
/*     read_nc({topo}, name, false); */
/*     fclose(fp); */
/*     backup_file(name); */
/*     fprintf(stdout, "%s .. ok\n", name); */
/*   } */

/* } */
