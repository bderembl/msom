%module pg
%include "common.i"
%include "grid/multigrid.i"
%include "predictor-corrector.i"
%include "poisson.i"
%include "timestep.i"

%apply (double * IN_ARRAY1, int DIM1) {
  (double * val1, int len1)
}
%apply (double * INPLACE_ARRAY1, int DIM1) {
  (double * val2, int len2)
}
%inline %{
  void pystep ( double * val1, int len1,
                double * val2, int len2);
%}
%apply (double * IN_ARRAY1, int DIM1) {
  (double * val3, int len3)
}



%inline %{
  void set_vars ();
  void pytrash_vars ();
  void pyset_contpar(int pycontpar);
  void pyadjust_contpar(double contpar_val);
  void pyinit_const(int pynl);
  void pyinit_last();
%}


