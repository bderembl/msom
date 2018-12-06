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
*  Routine originally written by A. Castillo, modified to get layered fields 
*/

struct InputMatrixl {
  // compulsory
  scalar * sl;
  FILE * fp;
  // optional
  int n;
  double ox, oy, width, oz;
  double arz;
};

void input_matrixl (struct InputMatrixl p) {
  for (scalar s in p.sl){
    if (p.width == 0.) p.width = L0;
    float width=0;
    fread(&width, sizeof(float), 1, p.fp);
    
    if (p.n != width) p.n = (int) width;
    
    float yp[p.n], xp[p.n], v[p.n][p.n];
    fread(&yp, sizeof(float), p.n, p.fp);
    for (int i = 0; i < p.n; i++) {
      fread(&xp[i], sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fread(&v[i][j], sizeof(float), 1, p.fp);
      }
    }
    foreach() {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < width)
        s[] = v[i][j];
      else
        s[] = 0.;
    }
  }
}

// modified output for scalar list
struct OutputMatrixl {
  scalar * fl;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_matrixl_nompi (struct OutputMatrixl p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n;
  float Delta = (float) L0/fn;
  for (scalar f in p.fl){
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
	v = interpolate (f, xp, yp);
      else {
	Point point = locate (xp, yp);
	assert (point.level >= 0);
	v = val(f);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  }
  fflush (p.fp);
}

void output_matrixl (struct OutputMatrixl p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));
  for (scalar f in p.fl){
  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
        field[i][j] = interpolate (f, xp, yp);
      }
      else {
        Point point = locate (xp, yp);
        field[i][j] = point.level >= 0 ? val(f) : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif


    fwrite (&fn, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < p.n; i++){
      float xp = Delta*i + X0 + Delta/2.;
      fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif
  }
  matrix_free (field);
}


/* /\* */
/* This function is an extension of output_matrix() inside "output.h" which allows */
/* to write 2D fields in a gnuplot-complatible format when running in MPI by */
/* performing a MPI_Reduce. */
/* *\/ */
/* void output_matrix_mpi (struct OutputMatrix p) */
/* { */
/*   if (p.n == 0) p.n = N; */
/*   if (!p.fp) p.fp = stdout; */
/*   float fn = p.n, Delta = L0/fn; */
/*   float ** field = matrix_new (p.n, p.n, sizeof(float)); */

/*   for (int i = 0; i < p.n; i++) { */
/*     float xp = Delta*i + X0 + Delta/2.; */
/*     for (int j = 0; j < p.n; j++) { */
/*       float yp = Delta*j + Y0 + Delta/2.; */
/*       if (p.linear) { */
/*         field[i][j] = interpolate (p.f, xp, yp); */
/*       } */
/*       else { */
/*         Point point = locate (xp, yp); */
/*         field[i][j] = point.level >= 0 ? val(p.f) : nodata; */
/*       } */
/*     } */
/*   } */

/*   if (pid() == 0) { // master */
/* @if _MPI */
/*     MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD); */
/* @endif */


/*     fwrite (&fn, sizeof(float), 1, p.fp); */
/*     for (int j = 0; j < p.n; j++) { */
/*       float yp = Delta*j + Y0 + Delta/2.; */
/*       fwrite (&yp, sizeof(float), 1, p.fp); */
/*     } */

/*     for (int i = 0; i < p.n; i++){ */
/*       float xp = Delta*i + X0 + Delta/2.; */
/*       fwrite (&xp, sizeof(float), 1, p.fp); */
/*       for (int j = 0; j < p.n; j++) { */
/*         fwrite (&field[i][j], sizeof(float), 1, p.fp); */
/*       } */
/*     } */
/*     fflush (p.fp); */
/*   } */
/* @if _MPI */
/*   else // slave */
/*   MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD); */
/* @endif */

/*   matrix_free (field); */
/* } */
