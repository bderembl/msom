
#ifdef _STOCHASTIC

int corrector_step = 0;
double tr_stoch = 0;
double itr_stoch = 0;
double amp_stoch = 1;
double nu_s = 0;

@define normal_noise() (sqrt(-2.*log(( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 2. ) ))*cos(2*pi*rand()/(double)RAND_MAX))

scalar var_s[];
scalar * n_stochl = NULL; // noise
scalar * s_stochl = NULL; // sigma of the noise


trace
double advection_pv(scalar * qol, scalar * qotl, scalar * pol, scalar * dqol, double dtmax)
{

  foreach() {
    double ju,jd;

    if (nl > 1){

      // upper layer
      int l = 0;
      scalar qo  = qol[l];
      scalar po  = pol[l];
      scalar pp  = ppl[l];
      scalar qot  = qotl[l];
#if _LS_RV
      scalar qp  = zetapl[l];
#endif
      scalar dqo = dqol[l];
      scalar po2  = pol[l+1];
      scalar pp2  = ppl[l+1];
      scalar s1 = strl[l];

      jd = jacobian(pp, po2) + jacobian(po, pp2);
      dqo[] += jacobian(pp, qo) + beta_effect(po) + s1[]*jd*idh1[l];
#if _LS_RV
      dqo[] += jacobian(po, qp);
#endif
      dqo[] += -qot[]*itr_stoch;

      // intermediate layers
      for (int l = 1; l < nl-1 ; l++) {
       
        qo  = qol[l];
        po  = pol[l];
        qot = qotl[l];
        pp  = ppl[l];
#if _LS_RV
        qp  = zetapl[l];
#endif
        dqo = dqol[l];
        po2  = pol[l+1];
        pp2  = ppl[l+1];
        scalar s0 = strl[l-1];
        scalar s1 = strl[l];

        ju = -jd;
        jd =  jacobian(pp, po2) + jacobian(po, pp2);
        dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju*idh0[l] + s1[]*jd*idh1[l];
#if _LS_RV
        dqo[] += jacobian(po, qp);
#endif
      dqo[] += -qot[]*itr_stoch;


      }

      // lower layer
      l = nl-1;

      qo  = qol[l];
      po  = pol[l];
      qot  = qotl[l];
      pp  = ppl[l];
#if _LS_RV
      qp  = zetapl[l];
#endif
      dqo = dqol[l];
      scalar s0 = strl[l-1];

      ju = -jd;
      dqo[] += jacobian(po, qo) + jacobian(pp, qo) + beta_effect(po) + s0[]*ju*idh0[l];
#if _LS_RV
      dqo[] += jacobian(po, qp);
#endif
      dqo[] += -qot[]*itr_stoch;

    }
    else{
      scalar dqo = dqol[0];
      dqo[] = 0.;
    }
  }

  // compute dtmax
   for (int l = 0; l < nl ; l++) {
     face vector uf[];
     scalar po = pol[l];
     comp_vel(po, uf);
     dtmax = timestep (uf, dtmax);
     po = ppl[l];
     comp_vel(po, uf);
     dtmax = timestep (uf, dtmax);
   }
  return dtmax;
}

/**
   generate normal noise with variance s_stoch^2
*/

double generate_noise(scalar * n_stochl, scalar *s_stochl)
{
  foreach() {
    for (int l = 0; l < nl; l++) {      
      scalar s_stoch  = s_stochl[l];
      scalar n_stoch  = n_stochl[l];
      n_stoch[] = amp_stoch*s_stoch[]*normal_noise();
    }
  }
}

static void advance_qg (scalar * output, scalar * input,
		        scalar * pol,
                        scalar * updates, double dt)
{
  
  corrector_step = (corrector_step+1)%2;
  float dts = sqrt(dt);
  if (corrector_step) {
    generate_noise(n_stochl, s_stochl);
    dts = dts/sqrt(2); // :to get sqrt(dt)/2 (in the predictor step, dt=dt/2)
    }
  
  double noiseV = 0., viscV = 0.;
  // biharmonic viscosity to tame the stochastic energy input
  foreach(reduction(+:noiseV) reduction(+:viscV)) { 
    noiseV += pol[]*nstochl[]*dh[];
    viscV  += pol[]*laplacian(laplacian(input[]))*dh[];
  }
  nu_s = noiseV/viscV;

  foreach() {
    for (int l = 0; l < nl ; l++) {
      scalar qi = input[l];
      scalar qo = output[l];
      scalar dq = updates[l];
      scalar n_stoch = n_stochl[l];
      qo[] = qi[] + dq[]*dt + n_stoch[]*dts + nu_s*laplacian(laplacian(qi[]));
    }
  }
  boundary(output);
}



void set_vars_stoch()
{
  fprintf(stdout,"Create stochastic variables .. ");

  int bc_type = 0;
  s_stochl   = create_layer_var(s_stochl,nl,bc_type);
  n_stochl   = create_layer_var(n_stochl,nl,bc_type);

  fprintf(stdout, ".. ok\n");

  fprintf(stdout, "Read stochastic sigma files (std deviation of the noise):\n");
  char name[80];
  FILE * fp;
  sprintf (name,"s_stoch_%dl_N%d.bas", nl,N);
  if ((fp = fopen (name, "r"))) {
    input_matrixl (s_stochl, fp);
    fclose(fp);
    fprintf(stdout, "%s .. ok\n", name);
  }
}


void trash_vars_stoch(){
  delete(s_stochl), free (s_stochl), s_stochl = NULL;
  delete(n_stochl), free (n_stochl), n_stochl = NULL;
}

event defaults (i = 0){
  set_vars_stoch();
}

event cleanup (i = end, last) {
  trash_vars_stoch();
}

#endif
