
%apply (double * IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * po_py, int len1, int len2, int len3)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * bfn_dp_py, int len4, int len5, int len6)
}

%inline %{
  void pystep_bfn ( double * po_py, int len1, int len2, int len3,
                double * bfn_dp_py, int len4, int len5, int len6,
                    double direction);
%}

%inline %{
  void set_vars_bfn ();
  void trash_vars_bfn ();
%}
