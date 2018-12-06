/**
   Planetary geostrophic model

   qcc -lm -O3 pg.c -o pg.e
*/

#include "grid/multigrid.h"
#include "../msqg/auxiliar_input.h"
#include "pg.h"

/**
 spatially varying non dimensional diffusivity coef. in dimensional units
$$
\kappa_v^* = \kappa \frac{N^2H^4}{\beta L^3}
$$
$$
\kappa_h^* = \kappa \frac{N^2H^4}{\beta L^3}*a^2*\frac{L^2}{H^2}
$$
 */
double k (double x, double y, double s) { return (3e-4);}

/**
   wind stress and wind stress derivative
in dimensional units [N/m^2]
$$
\tau_0^* = \tau \frac{\rho_0 N^2H^3}{2\pi L}
$$
*/

double tau0 = 0.12;

double taux   (double x, double y){ return (tau0*sin(2*(y-ys)*pi));}
double taux_y (double x, double y){ return (2*pi*tau0*cos(2*(y-ys)*pi));}
/* double taux   (double x, double y){ return (0.);} */
/* double taux_y (double x, double y){ return (0.);} */
double tauy   (double x, double y){ return (0.);}
double tauy_x (double x, double y){ return (0.);}




int main() {
  N = 64;
  nl = 30;

  CFL = 0.4;
  DT = 1.e-2;
  tend = 1000;
  dtout = 10;
  /**
     physical parameters : $a$ is the aspect ratio of the basin $H/L$.
     $r$ is the non dimensional friction coefficient. The dimensional
     friction is $r^* = r \beta L $. tau_surf is the surface
     relaxation time scale. */

  a = sqrt(3.0e-2/k(0,0,0)); 
  r = 0.02; 
  tau_surf = 3.0e-2;

  /**
     We use a pseudo SOR impletation because for small r, the
     ellicptic system is not diagonally dominant. For omega = 1, we
     recover the default basilisk implementation. omega<1 slows down
     the convergence rate, but at least it converges. */

  omega = 0.2;

  ys = 0.3;
  origin (0.0, ys);

  outdir = "outdir/";

  run(); 


}

event init (t = 0) {

  /**
     Initial conditions
  */
  FILE * fp ;
  if (fp = fopen ("b_init.bas", "r")){
    input_matrixl (bl, fp, oy=ys);
    fclose(fp);
  }

  if (fp = fopen ("u_init.bas", "r")) {
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

