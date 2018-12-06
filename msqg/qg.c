/**
# Multiple scale quasi geostrophic model

This is the driver file for qg.h.  We define the grid, the topography,
the initial conditions and the output routines.

compile with (openmp)
qcc -lm -llapacke qg.c -O3 -o qg.e -fopenmp
export OMP_NUM_THREADS=20 (?)
./qg.e

MPI:
CC99='mpicc -std=c99' qcc -D_MPI=1 -lm -O3 -llapacke qg.c -o qg.e
mpirun -np 16 ./qg.e
*/
#include "grid/multigrid.h"
#include "qg.h"
#include "auxiliar_input.h"


int main() {
  outdir = "outdir/";

/**
   Horizontal and vertical number of grid points */

  N = 512;

  char name[80];
  sprintf (name,"%sdh.bin", outdir);
  FILE * fp = fopen (name, "r");
  fseek(fp, 0, SEEK_END); 
  nl = ftell(fp)/sizeof(float); 
  fclose(fp);

  fprintf(stdout, "Config: N = %d, nl = %d\n",N, nl);

/**
   Time stepping parameters
*/
  DT = 5.e-2;
  CFL = 0.6;
  tend = 5000.0;
  dtout = 10.;
  dtfilter = 0.5;

/**
   Physical parameters: Size of the domain, Rossby number, Ekman
   number, Reynolds number, beta
*/
  Lt = 100; 
  Rom = 0.025;
  Ek = 1.0;
  Re = 15.0; //512:15 1024:50
  beta = 0.5;


  init_grid (N);
  size(Lt);
  run();
}


event init (i = 0) {


/**
   Layer thickness and large scale variables
*/
  char name[80];
  sprintf (name,"%sdh.bin", outdir);
  float dh[nl];
  FILE * fp = fopen (name, "r");  
  fread(&dh, sizeof(float), nl, fp);
  fclose(fp);

  for (int l = 0; l < nl ; l++)
    hl[l] = dh[l];

    
  sprintf (name,"%spsipg.bas0512", outdir);
  fp = fopen (name, "r");  
  input_matrixl (ppl, fp);
  fclose(fp);

  sprintf (name,"%sgppg.bas0512", outdir);
  fp = fopen (name, "r");  
  input_matrixl (gpl, fp);
  fclose(fp);
  
/**
   Initial conditions
*/
  foreach() 
    for (int l = 0; l < nl ; l++) {
      scalar qo  = qol[l];
      scalar po  = pol[l];
      qo[] = noise();
      po[] = 0.;
    }
  
}

event filter (t = 0; t <= tend+1e-10;  t += dtfilter) {
  fprintf(stdout,"Filter solution\n");
  invertq(pol,qol);
  wavelet_filter(pol, pofl, dtfilter);
  comp_q(pol,qol);
}

event writestdout (i++) {
/* event writestdout (i=1) { */
  scalar po = pol[0];
  scalar zeta = zetal[0];
  double ke = 0;
  foreach(reduction(+:ke))
    ke -= po[]*zeta[]*sq(Ro[]*Delta);

  fprintf (stdout,"i = %i, dt = %g, t = %g, ke = %g\n", i, dt, t, ke);
}

/**
   Write parameters
 */
event output (t = 0) {
  char name[80];
  sprintf (name,"%siBu.bas", outdir);
  FILE * fp = fopen (name, "w");
  output_matrixl (iBul, fp);
  fclose(fp);

  sprintf (name,"%ssig_filt.bas", outdir);
  fp = fopen (name, "w");
  output_matrixl ({sig_filt}, fp);
  fclose(fp);
}

event output (t = 0; t <= tend+1e-10;  t += dtout) {
  fprintf(stdout,"write file\n");

  if (t == 0) {
    // invert vorticity
    boundary(pol);
    invertq(pol,qol);
  }

  /**
   Rescale qo*/
  foreach()
    for (scalar qo in qol) 
      qo[] = qo[]*Ro[];
  boundary(qol);  

  char name[80];
  sprintf (name,"%spo%09d.bas", outdir, i);
  FILE * fp = fopen (name, "w");
  output_matrixl (pol, fp);
  fclose(fp);

  sprintf (name,"%sqo%09d.bas", outdir, i);
  fp = fopen (name, "w");
  output_matrixl (qol, fp);
  fclose(fp);

  sprintf (name,"%spf%09d.bas", outdir, i);
  fp = fopen (name, "w");
  output_matrixl (pofl, fp);
  fclose(fp);
  
  nbar = 0; // reset filter average

  /**
   Rescale qo*/
  foreach()
    for (scalar qo in qol)
      qo[] = qo[]/Ro[];
  boundary(qol);


  /* scalar l[]; */
  /* foreach() */
  /*   l[] = level; */
  /* sprintf (name,"%slevel%09d.dat", outdir, i); */
  /* fp = fopen (name, "w"); */
  /* output_field ({l}, fp); */
  /* fclose(fp); */


}

/* event adapt (t+=10) { */
/*  astats s = adapt_wavelet (pol, (double []){1e0, 1e0}, maxlevel = 9); */
/*  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc); */
/* } */
