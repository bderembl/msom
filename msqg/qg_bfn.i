
%apply (double * IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * varin_py, int len1, int len2, int len3)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * tend_py, int len4, int len5, int len6)
}

%inline %{
  void pystep_bfn ( double * varin_py, int len1, int len2, int len3,
                    double * tend_py, int len4, int len5, int len6,
                    double direction, int vartype);
%}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * po_py, int len7, int len8, int len9)
}

%apply (double * IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * qo_py, int len10, int len11, int len12)
}

%inline %{
  void pyq2p ( double * po_py, int len7, int len8, int len9,
               double * qo_py, int len10, int len11, int len12);
%}

%apply (double * IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * po_py, int len13, int len14, int len15)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * qo_py, int len16, int len17, int len18)
}

%inline %{
  void pyp2q ( double * po_py, int len13, int len14, int len15,
               double * qo_py, int len16, int len17, int len18);
%}


%inline %{
  void set_vars_bfn ();
  void trash_vars_bfn ();
%}
