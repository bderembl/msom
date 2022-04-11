/* input_matrix() reads a gnuplot-compatible binary file, i.e. single precision *
*  matrix stored in the following format:                                       *
*                                                                               *
*  <N+1>  <y0>   <y1>   <y2>  ...  <yN>                                         *
*  <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>                                         *
*  <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>                                         *
*                                                                               *
*  For example: input_matrix(T,fp,N,X0,Y0,L0);                                  *
*  Reads a square matrix of size N from file descriptor "fp" into scalar field  *
*  T, defined at points (X,Y).                                                  *
*                                                                               *
*  Routine originally written by A. Castillo, modified to get layered fields    *
*/

struct InputMatrixl {
  // compulsory
  scalar  s;
  FILE * fp;
  // optional
  int n;
  double ox, oy, width, oz;
  double arz;
};

void input_matrixl (struct InputMatrixl p) {
//  int il = 0;
  scalar s = p.s;
  foreach_layer()  {
    if (p.width == 0.) p.width = L0;
    float width=0;
    fread(&width, sizeof(float), 1, p.fp);
    
    if (p.n != width) p.n = (int) width;
    
    float yp[p.n], xp[p.n];
    float ** v = matrix_new (p.n, p.n, sizeof(float));
    fread(&yp, sizeof(float), p.n, p.fp);
    for (int i = 0; i < p.n; i++) {
      fread(&xp[i], sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fread(&v[i][j], sizeof(float), 1, p.fp);
      }
    }

    foreach(noauto) {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < width){
        s[] = v[i][j];
      }
      else
        s[] = 0.;
    }
  matrix_free (v);
  s.dirty = true;
  }
  boundary({s});
}
