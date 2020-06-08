/**
   Back and Forth Nudging routines
 */
scalar * bfn_tendl = NULL;
scalar * bfn_forcl = NULL;

void set_vars_bfn(){
  bfn_forcl = create_layer_var(bfn_forcl,nl,0);
  bfn_tendl = create_layer_var(bfn_tendl,nl,0);
}

void trash_vars_bfn(){
  delete(bfn_forcl), free(bfn_forcl), bfn_forcl = NULL;
  delete(bfn_tendl), free(bfn_tendl), bfn_tendl = NULL;
}


/**
Time step routine (stream function or PV)
*/
void pystep_bfn ( double * varin_py, int len1, int len2, int len3,
                  double * tend_py, int len4, int len5, int len6,
                  double direction, int vartype) {

  /**
     forward integration: direction = 1
     backward integration: direction = -1

     integration on p: vartype = 0
     integration on q: vartype = 1
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

  reset_layer_var(bfn_tendl);

  if (vartype == 0){
    fprintf (stdout,"temporary disabled psi tendency\n");

    /* pyset_field(pol,varin_py); */
    
    /* reset_layer_var(qol); // we put the tendency in qol */
    /* comp_del2(pol, zetal, 0., 1.0); */
    
    /* dtmax = advection(zetal, pol, qol, dtmax); */
    /* dissip(zetal, qol); */
    /* ekman_friction(zetal, qol); */
    
    /* // compute the stream function tendency */
    /* invertq(bfn_tendl, qol); */
    
    /* pyget_field(bfn_tendl, tend_py); */
  } else if (vartype == 1) {

    pyset_field(qol,varin_py);
    invertq(pol, qol);
    
    comp_del2(pol, zetal, 0., 1.0);
    dtmax = advection(zetal, qol, pol, bfn_tendl, dtmax);
    dissip(zetal, bfn_tendl);
    ekman_friction(zetal, bfn_tendl);
    if (flag_topo)
      bottom_topography(pol,bfn_tendl);

    pyget_field(bfn_tendl, tend_py);
  }
}

/**
Conversion routines
*/
void pyq2p ( double * po_py, int len7, int len8, int len9,
             double * qo_py, int len10, int len11, int len12) {

  reset_layer_var(pol);

  pyset_field(qol,qo_py);
  invertq(pol, qol);  
  pyget_field(pol, po_py);
}

void pyp2q ( double * po_py, int len13, int len14, int len15,
             double * qo_py, int len16, int len17, int len18) {

  reset_layer_var(qol);

  pyset_field(pol,po_py);
  comp_q(pol, qol);
  pyget_field(qol, qo_py);
}
