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
                    scalar * de_j1l, scalar * de_j2l, scalar * de_j3l, 
                    double dt, double ediag)
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
      scalar s1 = strl[l];
      
      jd_1 = jacobian(po, po2);
      jd_2 = jacobian(pp, po2);
      jd_3 = jacobian(po, pp2);
      jc   = jacobian(po, pp );

      de_j1[] += (jacobian(po, qo) + s1[]*jd_1*idh1[l]       )*(-po[]*dt*(1-ediag)+ediag);
      de_j2[] += (jacobian(pp, qo) + s1[]*(jd_2 + jc)*idh1[l])*(-po[]*dt*(1-ediag)+ediag);
      de_j3[] += (                   s1[]*(jd_3 - jc)*idh1[l])*(-po[]*dt*(1-ediag)+ediag);

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
        scalar s0 = strl[l-1];
        scalar s1 = strl[l];

        ju_1 = -jd_1;
        ju_2 = -jd_3; // swap
        ju_3 = -jd_2; // swap
        jd_1 = jacobian(po, po2);
        jd_2 = jacobian(pp, po2);
        jd_3 = jacobian(po, pp2);
        jc   = jacobian(po, pp );

        de_j1[] += (jacobian(po, qo) + s0[]*ju_1*idh0[l] + s1[]*jd_1*idh1[l]              )*(-po[]*dt*(1-ediag)+ediag);
        de_j2[] += (jacobian(pp, qo) + s0[]*(ju_2 + jc)*idh0[l] + s1[]*(jd_2 + jc)*idh1[l])*(-po[]*dt*(1-ediag)+ediag);
        de_j3[] += (                   s0[]*(ju_3 - jc)*idh0[l] + s1[]*(jd_3 - jc)*idh1[l])*(-po[]*dt*(1-ediag)+ediag);
      }

      // lower layer
      l = nl-1;

      qo  = qol[l];
      po  = pol[l];
      pp  = ppl[l];
      de_j1 = de_j1l[l];
      de_j2 = de_j2l[l];
      de_j3 = de_j3l[l];
      scalar s0 = strl[l-1];

      ju_1 = -jd_1;
      ju_2 = -jd_3; // swap
      ju_3 = -jd_2; // swap
      jc   = jacobian(po, pp );

      de_j1[] += (jacobian(po, qo) + s0[]*ju_1*idh0[l]       )*(-po[]*dt*(1-ediag)+ediag);
      de_j2[] += (jacobian(pp, qo) + s0[]*(ju_2 + jc)*idh0[l])*(-po[]*dt*(1-ediag)+ediag);
      de_j3[] += (                   s0[]*(ju_3 - jc)*idh0[l])*(-po[]*dt*(1-ediag)+ediag);
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
void dissip_de  (scalar * zetal, scalar * dqol, scalar * pol, double dt, double ediag)
{
  comp_del2(zetal, tmpl, 0., 1.);
  comp_stretch(zetal, tmp2l, 0., 1.);

  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar dqo = dqol[l];
      scalar p4 = tmpl[l];
      scalar str = tmp2l[l];
      scalar po = pol[l];
//      scalar pp  = ppl[l];
      dqo[] += (p4[] + str[])*iRe*(-po[]*dt*(1-ediag)+ediag);
      dqo[] += iRe4*laplacian(p4)
        *(-po[]*dt*(1-ediag)+ediag);
    }
  comp_stretch(tmpl, tmp2l, 0., 1.);
  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar dqo = dqol[l];
      scalar str = tmp2l[l];
      scalar po = pol[l];
//      scalar pp  = ppl[l];
      dqo[] += iRe4*(str[])
        *(-po[]*dt*(1-ediag)+ediag);
    }



}

trace
void ekman_friction_de (scalar * zetal, scalar * dqol, scalar * pol, double dt, double ediag)
{

  foreach(){
    scalar dqos = dqol[0];
    scalar zetas = zetal[0];
    scalar pos   = pol[0];

    scalar dqob = dqol[nl-1];
    scalar zetab = zetal[nl-1];
    scalar pob   = pol[nl-1];
    dqos[] -= Eks*zetas[]*(-pos[]*dt*(1-ediag)+ediag);
    dqob[] -= Ekb*zetab[]*(-pob[]*dt*(1-ediag)+ediag);
  }
}


trace
void filter_de (scalar * qol, scalar * pol, scalar * de_ftl, 
                scalar * po_mft, double dtflt, double ediag)
{
  // use tmp2 here because tmp is used in wavelet_filter
  int nbar = 0;
  wavelet_filter(qol, pol, tmp2l, -dtflt, nbar);
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar de_ft = de_ftl[l];
      scalar tmp = tmp2l[l];
//      scalar po = pol[l];
//      scalar pp   = ppl[l];
      scalar pm = po_mft[l];
      //de_ft[] += tmp[]*(-po[]-ediag*pp[])*Ro[]*Ro[];
      de_ft[] += tmp[]*(-pm[]*dtflt*(1-ediag)+ediag);
      pm[] = 0;
    }
  nme_ft = 0;
}

void energy_tend (scalar * pol, double dt)
{
  comp_del2(pol, zetal, 0., 1.0);
  advection_de(zetal, pol, de_j1l, de_j2l, de_j3l, dt, ediag);
  dissip_de(zetal, de_vdl, pol, dt, ediag);
  ekman_friction_de(zetal, de_bfl, pol, dt, ediag);

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
  delete(de_bfl), free(de_bfl), de_bfl = NULL;
  delete(de_vdl), free(de_vdl), de_vdl = NULL;
  delete(de_j1l), free(de_j1l), de_j1l = NULL;
  delete(de_j2l), free(de_j2l), de_j2l = NULL;
  delete(de_j3l), free(de_j3l), de_j3l = NULL;
  delete(de_ftl), free(de_ftl), de_ftl = NULL;
  delete(tmp2l),  free(tmp2l),  tmp2l = NULL;
  delete(po_mft), free(po_mft), po_mft = NULL;
}

/**
   Filter
*/
event filter (t = dtflt; t <= tend+1e-10;  t += dtflt) {
  if (ediag>-1)
    filter_de (qol, pol, de_ftl, po_mft, dtflt, ediag);
}


event defaults (i = 0){
 if (ediag>-1) {
   fprintf(stdout,"Create variables for energy diagnostics .. ");
   set_vars_energy();
   fprintf(stdout,"ok\n");
 }
}
event cleanup (i = end, last) {
 if (ediag>-1) trash_vars_energy();
}

event comp_diag (i++) {
 if (ediag>-1) energy_tend (pol, dt);
}

/**
   Python interface routines (should be in .i file but I need foreach)
 */

void pystep_de ( double * po_py, int len1, int len2, int len3,
              double * de_bf_py, int len4, int len5, int len6,
              double * de_vd_py, int len7, int len8, int len9,
              double * de_j1_py, int len10, int len11, int len12,
              double * de_j2_py, int len13, int len14, int len15,
              double * de_j3_py, int len16, int len17, int len18,
              double * de_ft_py, int len19, int len20, int len21,
              int onlyKE) {

  int ediag = 1;
  double dt = 1.;

  pyset_field(pol,po_py);

  reset_layer_var(de_bfl);
  reset_layer_var(de_vdl);
  reset_layer_var(de_j1l);
  reset_layer_var(de_j2l);
  reset_layer_var(de_j3l);
  reset_layer_var(de_ftl);

  comp_del2(pol, zetal, 0., 1.0);
  comp_q(pol,qol);

  if (onlyKE == 1) {
  foreach()
    for (int l = 0; l < nl-1 ; l++) {
      scalar s = strl[l];
      s[] = 0.;
    }
  }

  advection_de(zetal, pol, de_j1l, de_j2l, de_j3l, dt, ediag);
  dissip_de(zetal, de_vdl, pol, dt, ediag);
  ekman_friction_de(zetal, de_bfl, pol, dt, ediag);
  //todo: wrong filter if onlyKE = 1
  filter_de (qol, pol, de_ftl, pol, dtflt, ediag);
  

  pyget_field(de_bfl,de_bf_py);
  pyget_field(de_vdl,de_vd_py);
  pyget_field(de_j1l,de_j1_py);
  pyget_field(de_j2l,de_j2_py);
  pyget_field(de_j3l,de_j3_py);
  pyget_field(de_ftl,de_ft_py);

}
