
%apply (double * IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * qo_py, int len1, int len2, int len3)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * de_bf_py, int len4, int len5, int len6)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * de_vd_py, int len7, int len8, int len9)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * de_j1_py, int len10, int len11, int len12)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * de_j2_py, int len13, int len14, int len15)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * de_j3_py, int len16, int len17, int len18)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * de_ft_py, int len19, int len20, int len21)
}

%inline %{
  void pystep ( double * qo_py, int len1, int len2, int len3,
                double * de_bf_py, int len4, int len5, int len6,
                double * de_vd_py, int len7, int len8, int len9,
                double * de_j1_py, int len10, int len11, int len12,
                double * de_j2_py, int len13, int len14, int len15,
                double * de_j3_py, int len16, int len17, int len18,
                double * de_ft_py, int len19, int len20, int len21
                );
%}

%inline %{
  void set_vars_energy ();
  void trash_vars_energy ();
%}
