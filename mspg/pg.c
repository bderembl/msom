/**
   Planetary geostrophic model

   qcc -lm -O3 pg.c -o pg.e
*/

#include "grid/multigrid.h"
#include "../msqg/auxiliar_input.h"
#include "pg.h"

// for mkdir
#include <sys/stat.h>
#include <sys/types.h>


/**
 spatially varying non dimensional diffusivity coef. in dimensional units
$$
\kappa_v^* = \kappa \frac{N^2H^4}{\beta L^3}
$$
$$
\kappa_h^* = \kappa \frac{N^2H^4}{\beta L^3}*a^2*\frac{L^2}{H^2}
$$
 */
 //shape function, multipy by kd
double k (double x, double y, double s) { return (1.0);}

/**
   wind stress and wind stress derivative
in dimensional units [N/m^2]
$$
\tau_0^* = \tau \frac{\rho_0 N^2H^3}{2\pi L}
$$
*/

/* double taux   (double x, double y){ return (sin(2*(y-ys)*pi));} */
double taux   (double x, double y){ return (0.);}
double taux_y (double x, double y){ return (2*pi*y*cos(2*(y-ys)*pi));}
/* double taux_y (double x, double y){ return (0.);} */
double tauy   (double x, double y){ return (0.);}
double tauy_x (double x, double y){ return (0.);}


int N0 = 64;


int main() {

/**
   Read input parameters
 */

  FILE * fp;
  if (fp = fopen("params.in", "rt")) {
    char tempbuff[100];
    char tmps1[80];
    char tmps2[80];
    char tmps3[80];

    while(fgets(tempbuff,100,fp)) {
      sscanf(tempbuff, "%15s = %15s # %15s", tmps1, tmps2, tmps3);
      if      (strcmp(tmps1,"N")    ==0) { N0    = atoi(tmps2); }
      else if (strcmp(tmps1,"nl")   ==0) { nl    = atoi(tmps2); }
      else if (strcmp(tmps1,"a")    ==0) { a     = atof(tmps2); }
      else if (strcmp(tmps1,"r")    ==0) { r     = atof(tmps2); }
      else if (strcmp(tmps1,"kd")   ==0) { kd    = atof(tmps2); }
      else if (strcmp(tmps1,"tau_s")==0) { tau_s = atof(tmps2); }
      else if (strcmp(tmps1,"tau0") ==0) { tau0  = atof(tmps2); }
      else if (strcmp(tmps1,"ys")   ==0) { ys    = atof(tmps2); }
      else if (strcmp(tmps1,"omega")==0) { omega = atof(tmps2); }
      else if (strcmp(tmps1,"DT")   ==0) { DT    = atof(tmps2); }
      else if (strcmp(tmps1,"CFL")  ==0) { CFL   = atof(tmps2); }
      else if (strcmp(tmps1,"tend") ==0) { tend  = atof(tmps2); }
      else if (strcmp(tmps1,"dtout")==0) { dtout = atof(tmps2); }
    }
    fclose(fp);
  } else {
    fprintf(stdout, "file params.in not found\n");
    exit(0);
  }

  N = N0;

  /**
     physical parameters : $a$ is the aspect ratio of the basin $H/L$.
     $r$ is the non dimensional friction coefficient. The dimensional
     friction is $r^* = r \beta L $. tau_surf is the surface
     relaxation time scale. */

  /**
     We use a pseudo SOR impletation because for small r, the
     ellicptic system is not diagonally dominant. For omega = 1, we
     recover the default basilisk implementation. omega<1 slows down
     the convergence rate, but at least it converges. */

  origin (0.0, ys);

  /**
     Create output directory and copy input parameter file for backup
  */

  char ch;
  char name[80];
  int idir;
  if (pid() == 0) { // master
    for (idir=1; idir<10000; idir++) {
      sprintf(outdir, "outdir_%04d/", idir);
      if (mkdir(outdir, 0777) == 0) {
        fprintf(stdout,"Writing output in %s\n",outdir);
        break;
      }
    }
  }
@if _MPI 
  MPI_Bcast(&idir, 1, MPI_INT, 0, MPI_COMM_WORLD);
 sprintf(outdir, "outdir_%04d/", idir);
@endif

  sprintf (name,"%sparams.in", outdir);
  FILE * source = fopen("params.in", "r");
  FILE * target = fopen(name, "w");
  while ((ch = fgetc(source)) != EOF)
    fputc(ch, target);
  fclose(source);
  fclose(target);
  
  run(); 


}

event init (t = 0) {

  /**
     Initial conditions
  */
  FILE * fp ;
  if (fp = fopen ("b0.bas", "r")){
    input_matrixl (bl, fp, oy=ys);
    fclose(fp);
  }

  if (fp = fopen ("u0.bas", "r")) {
  input_matrixl ((scalar *) ul, fp, oy=ys);
  fclose(fp);
  }
  /**
     QG forcing
  */
  if (fp = fopen ("bf_pg.bas", "r")) {
  input_matrixl (b_forcl, fp, oy=ys);
  fclose(fp);
  }

  /**
     Surface forcing
  */
  foreach()
    b_surf[] = 6*cos(pi*(y-ys));

}

event writestdout (i++) {
  printf ("i = %i, dt = %g, t = %g\n", i, dt, t);
}


/* event adapt (i++) { */
/*  astats s = adapt_wavelet (bl, (double []){1.0}, maxlevel = 10); */
/*  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc); */
/* } */


/**
   write output: first compute running average
 */

event comp_output(i+=10) {
  foreach()
    for (int l = 1; l < nl+1; l++) {
      scalar b = bl[l];
      scalar bm = b_mel[l];
      bm[] = (bm[]*nme + b[])/(nme + 1);
    }
  boundary (b_mel);


  foreach_face(){
    for (int l = 1; l < nl+1 ; l++) {
      face vector u = ul[l];
      face vector um = u_mel[l];
      um.x[] = (um.x[]*nme + u.x[])/(nme + 1);
    }
  }

  for (int l = 1; l < nl+1 ; l++) {
    face vector um = u_mel[l];
    boundary ((scalar *){um});
  }
  nme++;
}
event writeconst (t = 0) {
  char name[80];
  sprintf (name,"%spsibt.bas", outdir, i);
  FILE * fp = fopen (name, "w");
  output_matrixl ({psibt}, fp);
  fclose(fp);
}


event writestate (t = 0; t <= tend+1e-10;  t += dtout) {
  printf ("i = %d t = %g\n", i, t);

  char name[80];
  sprintf (name,"%sb%09d.bas", outdir, i);
  FILE * fp = fopen (name, "w");
  output_matrixl (b_mel, fp);
  fclose(fp);

  sprintf (name,"%su%09d.bas", outdir, i);
  fp = fopen (name, "w");
  output_matrixl ((scalar *) u_mel, fp);
  fclose(fp);

  nme = 0; // reset running average
}

