/**
   Energy dynagnostics of the MSQG model

   We multiply all terms of the PV equation by -po*dt
*/

scalar * de_bfl = NULL;
scalar * de_vdl = NULL;
scalar * de_j1l = NULL;
scalar * de_j2l = NULL;
scalar * de_j3l = NULL;
scalar * de_ftl = NULL;
scalar * tmp2l = NULL;

scalar * po_mft = NULL;

int nme_ft = 0;

/**
   J1 = j(psi, q)
   J2 = j(psi_pg, q)
   J3 = j(psi, q_pg)
 */
trace
void advection_de  (scalar * qol, scalar * pol,
                    scalar * de_j1l, scalar * de_j2l, scalar * de_j3l, double dt)
{
  foreach() {
    double ju_1,jd_1;
    double ju_2,jd_2;
    double ju_3,jd_3;
    double jc;
    
    if (nl > 1){

      // upper layer
      int l = 0;
      scalar qo  = qol[l];
      scalar po  = pol[l];
      scalar pp  = ppl[l];
      scalar po2 = pol[l+1];
      scalar pp2 = ppl[l+1];
      scalar de_j1 = de_j1l[l];
      scalar de_j2 = de_j2l[l];
      scalar de_j3 = de_j3l[l];
      scalar s1 = str1l[l];
      
      jd_1 = jacobian(po, po2);
      jd_2 = jacobian(pp, po2);
      jd_3 = jacobian(po, pp2);
      jc   = jacobian(po, pp );

      de_j1[] += (jacobian(po, qo) + s1[]*jd_1          )*(-po[]-ediag*pp[])*dt;
      de_j2[] += (jacobian(pp, qo) + s1[]*jd_2 + s1[]*jc)*(-po[]-ediag*pp[])*dt;
      de_j3[] += (                   s1[]*jd_3 - s1[]*jc)*(-po[]-ediag*pp[])*dt;

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        qo  = qol[l];
        po  = pol[l];
        pp  = ppl[l];
        po2 = pol[l+1];
        pp2 = ppl[l+1];
        de_j1 = de_j1l[l];
        de_j2 = de_j2l[l];
        de_j3 = de_j3l[l];
        scalar s0 = str0l[l];
        scalar s1 = str1l[l];

        ju_1 = -jd_1;
        ju_2 = -jd_3; // swap
        ju_3 = -jd_2; // swap
        jd_1 = jacobian(po, po2);
        jd_2 = jacobian(pp, po2);
        jd_3 = jacobian(po, pp2);
        jc   = jacobian(po, pp );

        de_j1[] += (jacobian(po, qo) + s0[]*ju_1 + s1[]*jd_1                 )*(-po[]-ediag*pp[])*dt;
        de_j2[] += (jacobian(pp, qo) + s0[]*ju_2 + s1[]*jd_2 + (s1[]+s0[])*jc)*(-po[]-ediag*pp[])*dt;
        de_j3[] += (                   s0[]*ju_3 + s1[]*jd_3 - (s1[]+s0[])*jc)*(-po[]-ediag*pp[])*dt;
      }

      // lower layer
      l = nl-1;

      qo  = qol[l];
      po  = pol[l];
      pp  = ppl[l];
      de_j1 = de_j1l[l];
      de_j2 = de_j2l[l];
      de_j3 = de_j3l[l];
      scalar s0 = str0l[l];

      ju_1 = -jd_1;
      ju_2 = -jd_3; // swap
      ju_3 = -jd_2; // swap
      jc   = jacobian(po, pp );

      de_j1[] += (jacobian(po, qo) + s0[]*ju_1          )*(-po[]-ediag*pp[])*dt;
      de_j2[] += (jacobian(pp, qo) + s0[]*ju_2 + s0[]*jc)*(-po[]-ediag*pp[])*dt;
      de_j3[] += (                   s0[]*ju_3 - s0[]*jc)*(-po[]-ediag*pp[])*dt;
    }
    else{
      scalar de_j1 = de_j1l[0];
      scalar de_j2 = de_j2l[0];
      scalar de_j3 = de_j3l[0];
      de_j1[] = 0;
      de_j2[] = 0;
      de_j3[] = 0;
    }
  }
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
      dqo[] += p4[]*iRe*(-po[]-ediag*pp[])*dt;
      dqo[] += iRe4*(p4[1] + p4[-1] + p4[0,1] + p4[0,-1] - 4*p4[])/(sq(Delta))
        *(-po[]-ediag*pp[])*dt;
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
    dqo[] -= Ek*zeta[]*(-po[]-ediag*pp[])*dt;
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
      de_ft[] += tmp[]*(-pm[]-ediag*pp[]);
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
}

void set_vars_energy(){
  de_bfl = create_layer_var(de_bfl,nl,0);
  de_vdl = create_layer_var(de_vdl,nl,0);
  de_j1l = create_layer_var(de_j1l,nl,0);
  de_j2l = create_layer_var(de_j2l,nl,0);
  de_j3l = create_layer_var(de_j3l,nl,0);
  de_ftl = create_layer_var(de_ftl,nl,0);
  tmp2l  = create_layer_var(tmp2l, nl,0);
  po_mft = create_layer_var(po_mft, nl,0);
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
