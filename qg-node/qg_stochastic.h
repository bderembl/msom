#ifdef _STOCHASTIC

double amp_stoch = 0;
double L_filt = 0;
int corrector_step = 0;

// random noise define on scalar instead of vertex in order to use the wavelet transform

scalar sig_lev[];
scalar n_stoch[];

@define normal_noise() (sqrt(-2.*log(( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 2. ) ))*cos(2*pi*rand()/(double)RAND_MAX))

event init_stoch (i = 0) {

  reset ({sig_lev}, 0.);

  // low pass filter
  for (int l = depth(); l >= 0; l--) {
    foreach_level (l) {
      double ref_flag = 0;
      if (l < depth())
        foreach_child()
          ref_flag += sig_lev[];
      if (ref_flag > 0)
        sig_lev[] = 1;
      else{

        if (L_filt > 2*Delta)
          sig_lev[] = 0;
        else if (L_filt <= 2*Delta && L_filt > Delta)
          sig_lev[] = 1-(L_filt-Delta)/Delta;
        else
          sig_lev[] = 1;
      }
    }
    boundary_level ({sig_lev}, l);
  }

  // high pass filter
  for (int l = depth(); l >= 0; l--) {
    foreach_level (l)
      sig_lev[] = 1 - sig_lev[];
    boundary_level ({sig_lev}, l);
  }

}

void generate_noise(scalar psi, double amp_stoch)
{
  foreach() {
    psi[] = amp_stoch*normal_noise();
  }
  scalar w[];
  
  wavelet(psi,w);
  for (int l = 0; l <= depth(); l++) {
    foreach_level (l)
      w[] *= sig_lev[];
    boundary_level ({w}, l);
  }

  inverse_wavelet (psi, w);
  
}

// endif stochastic
#endif 
