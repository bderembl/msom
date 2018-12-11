%{
  struct InputMatrixl {
    // compulsory
    scalar * sl;
    FILE * fp;
    // optional
    int n;
    double ox, oy, width, oz;
    double arz;
  };
  
  extern void input_matrixl (struct InputMatrixl p);
%}
