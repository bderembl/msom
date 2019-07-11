/**
   Energy dynagnostics of the MSQG model

   We multiply all terms of the PV equation by po/f
   We also multiply the entire equation by 1/f
*/

scalar * de_bfl = NULL;
scalar * de_vdl = NULL;
scalar * de_j1l = NULL;
scalar * de_j2l = NULL;
scalar * de_j3l = NULL;
scalar * de_ftl = NULL;
scalar * tmp2l = NULL;

scalar * po_mft = NULL;

// temporary
scalar * qol_prev = NULL;
scalar * de_tol = NULL;

int nme_ft = 0;

trace
void combine_jac_de(scalar * jac1l, scalar * jacal, double add,
                 double p0, double p1, scalar * pol)
{
  foreach() {
    if (nl > 1){
      // upper layer
      int l = 0;
      scalar jac1 = jac1l[l];
      scalar jaca = jacal[l];
      scalar Fr1 = Frl[l];
      scalar po_r  = pol[l];
      scalar pp_r  = ppl[l];
      double b1 = sq(Fr1[]/Ro[])/( dhc[l]*dhf[l]);

      jaca[] = add*jaca[] + p1*b1*jac1[]*(-po_r[]-ediag*pp_r[])*dt*sq(Ro[]);

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        scalar jac0 = jac1l[l-1];
        scalar jac1 = jac1l[l];
        scalar jaca = jacal[l];
        scalar po_r1 = pol[l];
        scalar pp_r1 = ppl[l];
        
        scalar Fr0 = Frl[l-1];
        scalar Fr1 = Frl[l];
        
        double b0 = sq(Fr0[]/Ro[])/( dhc[l-1]*dhf[l]);
        double b1 = sq(Fr1[]/Ro[])/( dhc[l]*dhf[l]);

        jaca[] = add*jaca[] + (- p0*b0*jac0[] + p1*b1*jac1[])
          *(-po_r1[]-ediag*pp_r1[])*dt*sq(Ro[]);

      }
      // lower layer
      l = nl-1;
      jac1 = jac1l[l-1];
      jaca = jacal[l];
      po_r  = pol[l];
      pp_r  = ppl[l];
      scalar Fr0 = Frl[l-1];
      double b0 = sq(Fr0[]/Ro[])/( dhc[l-1]*dhf[l]);

      jaca[] = add*jaca[] - p0*b0*jac1[]*(-po_r[]-ediag*pp_r[])*dt*sq(Ro[]);

    }
    else{
      scalar jaca = jacal[0];
      jaca[] = 0.;
    }
  }
}

trace
void comp_part_stretch_de(scalar * pol, scalar * stretchl, double add, double fac)
{

  foreach() {

    if (nl > 1){
      // upper layer
      int l = 0;
      scalar po_1 = pol[l];
      scalar po_2 = pol[l+1];
      scalar stretch = stretchl[l];
      scalar Fr1 = Frl[l];
      double b1 = sq(Fr1[]/Ro[])/( dhc[l]*dhf[l]);

      stretch[] = add*stretch[] - fac*b1*po_1[] ;

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        scalar po_0 = pol[l-1];
        scalar po_1 = pol[l];
        scalar po_2 = pol[l+1];
        scalar stretch = stretchl[l];
        
        scalar Fr0 = Frl[l-1];
        scalar Fr1 = Frl[l];
        
        double b0 = sq(Fr0[]/Ro[])/( dhc[l-1]*dhf[l]);
        double b1 = sq(Fr1[]/Ro[])/( dhc[l]*dhf[l]);

        stretch[] = add*stretch[] - fac*(b0 + b1)*po_1[] ;

      }
      // lower layer
      l = nl-1;
      scalar po_0 = pol[l-1];
      po_1 = pol[l];
      stretch = stretchl[l];
      scalar Fr0 = Frl[l-1];
      double b0 = sq(Fr0[]/Ro[])/( dhc[l-1]*dhf[l]);

      stretch[] = add*stretch[] - fac*b0*po_1[] ;

    }
    else{
      scalar stretch = stretchl[0];
      stretch[] = 0.;
    }
  }

  boundary(stretchl);
}







trace
void jacobian_de(scalar po, scalar qo, scalar jac, double add,
                 scalar po_r, scalar pp_r, double dt)
{
  foreach()
    jac[] = add*jac[] +
    (( qo[1, 0 ]-qo[-1, 0])*(po[0, 1]-po[ 0 ,-1])
     +(qo[0 ,-1]-qo[ 0 ,1])*(po[1, 0]-po[-1, 0 ])
     + qo[1, 0 ]*( po[1,1 ] - po[1,-1 ])
     - qo[-1, 0]*( po[-1,1] - po[-1,-1])
     - qo[ 0 ,1]*( po[1,1 ] - po[-1,1 ])
     + qo[0 ,-1]*( po[1,-1] - po[-1,-1])
     + po[ 0 ,1]*( qo[1,1 ] - qo[-1,1 ])
     - po[0 ,-1]*( qo[1,-1] - qo[-1,-1])
     - po[1, 0 ]*( qo[1,1 ] - qo[1,-1 ])
     + po[-1, 0]*( qo[-1,1] - qo[-1,-1]))
    *Ro[]/(12.*Delta*Delta)
    *(-po_r[]-ediag*pp_r[])*dt*sq(Ro[]);
}

/**
   J1 = j(psi, q)
   J2 = j(psi_pg, q)
   J3 = j(psi, q_pg)
 */
trace
void advection_de  (scalar * qol, scalar * pol, 
                    scalar * de_j1l, scalar * de_j2l, scalar * de_j3l, double dt)
{  
  for (int l = 0; l < nl ; l++) {
    scalar qo  = qol[l];
    scalar po  = pol[l];
    scalar pp  = ppl[l];
    scalar de_j1 = de_j1l[l];
    scalar de_j2 = de_j2l[l];
    jacobian_de(po, qo, de_j1, 1., po, pp, dt); // -J(p_qg, zeta)
    jacobian_de(pp, qo, de_j2, 1., po, pp, dt); // -J(p_pg, zeta)
    /* jacobian_de(po, qo, de_j1, 1.); // -J(p_qg, zeta) */
    /* jacobian_de(pp, qo, de_j2, 1.); // -J(p_pg, zeta) */
  }
  
   for (int l = 0; l < nl-1 ; l++) {
    scalar po1  = pol[l];
    scalar po2  = pol[l+1];
    scalar pp  = ppl[l];
    scalar jac1 = tmpl[l];
    /* jacobian_de(po1, po2, jac1, 0., po1, pp, dt); // -J(p_l, p_l+1) */
    jacobian(po1, po2, jac1, 0.); // -J(p_l, p_l+1)
  }
 combine_jac_de(tmpl, de_j1l, 1., 1. , 1., pol);

   for (int l = 0; l < nl-1 ; l++) {
    scalar po1  = pol[l];
    scalar pp1  = ppl[l];
    scalar po2  = pol[l+1];
    scalar jac1 = tmpl[l];
    /* jacobian_de(pp1, po2, jac1, 0., po1, pp1, dt); // -J(pp_l, p_l+1) */
    jacobian(pp1, po2, jac1, 0.); // -J(pp_l, p_l+1)
  }
 combine_jac_de(tmpl, de_j2l, 1., 0., 1., pol);
 combine_jac_de(tmpl, de_j3l, 1., 1., 0., pol);


   for (int l = 0; l < nl-1 ; l++) {
    scalar po1  = pol[l];
    scalar jac1 = tmpl[l];
    scalar pp1  = ppl[l];
    scalar pp2  = ppl[l+1];
    /* jacobian_de(po1, pp2, jac1, 0., po1, pp1, dt); // -J(p_l, pp_l+1) */
    jacobian(po1, pp2, jac1, 0.); // -J(p_l, pp_l+1)
  }
 combine_jac_de(tmpl, de_j2l, 1., 1., 0., pol);
 combine_jac_de(tmpl, de_j3l, 1., 0., 1., pol);

   for (int l = 0; l < nl ; l++) {
    scalar po1  = pol[l];
    scalar pp1  = ppl[l];
    scalar jac1 = tmpl[l];
    /* jacobian_de(pp1, po2, jac1, 0., po1, pp1, dt); // -J(pp_l, p_l+1) */
    jacobian_de(po1, pp1, jac1, 0., po1, pp1, dt); // -J(p_l, pp_l)
  }
   // add and subtract J(p,pp) in j2 and j3
   comp_part_stretch_de(tmpl, de_j3l, 1.,  1.);
   comp_part_stretch_de(tmpl, de_j2l, 1., -1.);

}


trace
void dissip_de  (scalar * zetal, scalar * dqol, scalar * pol, double dt)
{
  double iRe = 1/Re;
  double iRe4 = -1/Re4;
  comp_del2(zetal, tmpl, 0., 1.);

  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar dqo = dqol[l];
      scalar p4 = tmpl[l];
      scalar po = pol[l];
      scalar pp  = ppl[l];
      dqo[] += p4[]*iRe*(-po[]-ediag*pp[])*dt*sq(Ro[]);
      dqo[] += iRe4*(p4[1] + p4[-1] + p4[0,1] + p4[0,-1] - 4*p4[])/(sq(Delta))
        *(-po[]-ediag*pp[])*dt*sq(Ro[]);
    }
}

trace
void bottom_friction_de (scalar * zetal, scalar * dqol, scalar * pol, double dt)
{
  scalar dqo  = dqol[nl-1];
  scalar zeta = zetal[nl-1];
  scalar po   = pol[nl-1];
  scalar pp   = ppl[nl-1];

  foreach()
    dqo[] -= Ek*zeta[]*(-po[]-ediag*pp[])*dt*sq(Ro[]);
}


trace
void filter_de (scalar * qol, scalar * pol, scalar * de_ftl)
{
  // use tmp2 here because tmp is used in wavelet_filter
  int nbar = 0;
  wavelet_filter(qol, pol, tmp2l, -1.0, nbar);
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar de_ft = de_ftl[l];
      scalar tmp = tmp2l[l];
      scalar po = pol[l];
      scalar pp   = ppl[l];
      scalar pm = po_mft[l];
      // no mutliplication by dt here
      //de_ft[] += tmp[]*(-po[]-ediag*pp[])*Ro[]*Ro[];
      de_ft[] += tmp[]*(-pm[]-ediag*pp[])*sq(Ro[]);
      pm[] = 0;
    }
  nme_ft = 0;
}

void energy_tend (scalar * pol, double dt)
{
  comp_del2(pol, zetal, 0., 1.0);
  advection_de(zetal, pol, de_j1l, de_j2l, de_j3l, dt);
  dissip_de(zetal, de_vdl, pol, dt);
  bottom_friction_de(zetal, de_bfl, pol, dt);
  // filter part in qg.c

  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar po = pol[l];
      scalar pm = po_mft[l];
      pm[] = (pm[]*nme_ft + po[])/(nme_ft+1);
    }
  nme_ft += 1;


  // temporary
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar qo1 = qol[l];
      scalar qo0 = qol_prev[l];
      scalar po = pol[l];
      scalar pp   = ppl[l];
      scalar de_to = de_tol[l];
      de_to[] += (qo1[] - qo0[])/dt*(-po[]-ediag*pp[])*dt*sq(Ro[]);
      qo0[] = qo1[];
    }

}

void set_vars_energy(){
  de_bfl = create_layer_var(de_bfl,nl);
  de_vdl = create_layer_var(de_vdl,nl);
  de_j1l = create_layer_var(de_j1l,nl);
  de_j2l = create_layer_var(de_j2l,nl);
  de_j3l = create_layer_var(de_j3l,nl);
  de_ftl = create_layer_var(de_ftl,nl);
  tmp2l  = create_layer_var(tmp2l, nl);
  po_mft = create_layer_var(po_mft, nl);
  // temporary
  qol_prev  = create_layer_var(qol_prev, nl);
  de_tol = create_layer_var(de_tol,nl);

}

void trash_vars_energy(){
  free(de_bfl), de_bfl = NULL;
  free(de_vdl), de_vdl = NULL;
  free(de_j1l), de_j1l = NULL;
  free(de_j2l), de_j2l = NULL;
  free(de_j3l), de_j3l = NULL;
  free(de_ftl), de_ftl = NULL;
  free(tmp2l), tmp2l = NULL;
  free(po_mft), po_mft = NULL;
  // temporary
  free(qol_prev), qol_prev = NULL;
  free(de_tol), de_tol = NULL;
}

event defaults (i = 0){
 if (ediag>-1) set_vars_energy();
}
event cleanup (i = end, last) {
 if (ediag>-1) trash_vars_energy();
}

event comp_diag (i++) {
 if (ediag>-1) energy_tend (pol, dt);
}
