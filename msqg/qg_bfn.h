/**
   Back and Forth Nudging routines
 */
scalar * bfn_dpl = NULL;
scalar * bfn_forcl = NULL;

void set_vars_bfn(){
  bfn_forcl = create_layer_var(bfn_forcl,nl,0);
  bfn_dpl = create_layer_var(bfn_dpl,nl,0);
}

void trash_vars_bfn(){
  delete(bfn_forcl), free(bfn_forcl), bfn_forcl = NULL;
  delete(bfn_dpl), free(bfn_dpl), bfn_dpl = NULL;
}


void pystep_bfn ( double * po_py, int len1, int len2, int len3,
              double * bfn_dp_py, int len4, int len5, int len6,
                  double direction) {

  /**
     forward integration: direction = 1
     backward integration: direction = -1
   */
  double dtmax = DT;

  if (direction > 0) {
    if (Re  == 0) iRe  = 0.; else iRe  =  1/Re;
    if (Re4 == 0) iRe4 = 0.; else iRe4 = -1/Re4;
    Eks = fabs(Eks);
    Ekb = fabs(Ekb);
  } else {
    if (Re  == 0) iRe  = 0.; else iRe  =  -1/Re;
    if (Re4 == 0) iRe4 = 0.; else iRe4 = 1/Re4;
    Eks = -fabs(Eks);
    Ekb = -fabs(Ekb);
  }

  pyset_field(pol,po_py);

  reset_layer_var(bfn_dpl);
  reset_layer_var(qol); // we put the tendency in qol 

  comp_del2(pol, zetal, 0., 1.0);

  dtmax = advection(zetal, pol, qol, dtmax);
  dissip(zetal, qol);
  ekman_friction(zetal, qol);

  // compute the stream function tendency
  invertq(bfn_dpl, qol);
  
  pyget_field(bfn_dpl, bfn_dp_py);

}
