
%apply (double * IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * val1, int len1, int len2, int len3)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * val2, int len4, int len5, int len6)
}

%inline %{
  void pystep ( double * val1, int len1, int len2, int len3,
                double * val2, int len4, int len5, int len6,
                char * id);
%}
