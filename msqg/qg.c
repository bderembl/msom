/**
# Multiple scale quasi geostrophic model

This is the driver file for qg.h.  We define the grid, the topography,
the initial conditions and the output routines.

compile with (openmp)
qcc -lm qg.c -O3 -o qg.e -fopenmp (-llapacke)
export OMP_NUM_THREADS=20 (?)
./qg.e

MPI:
CC99='mpicc -std=c99' qcc -D_MPI=1 -lm -O3 qg.c -o qg.e -grid=multigrid (-llapacke)
mpirun -np 16 ./qg.e

HPC:
qcc -D_MPI=1 -grid=multigrid -source qg.c     (-D MKL)
rsync _qg.c
mpicc -Wall -std=c99 -O2 _qg.c -lm -o qg.e    (-mkl)
*/

#include "grid/multigrid.h"
#include "auxiliar_input.h"
#include "qg.h"
#include "qg_energy.h"
#include "qg_bfn.h"

int main(int argc,char* argv[]) {

  // Search for the configuration file with a given path or read params.in 
  if (argc == 2) {
    read_params(argv[1]);
  } else {
    read_params("params.in");
  }

  create_outdir();

  init_grid (N);
  size(L0);
  run();
}

/**
   Initial conditions
*/
event init (i = 0) {

  FILE * fp;
  if ((fp = fopen("p0.bas", "r"))) {
    input_matrixl (pol, fp);
    fclose(fp);
  } else {
  foreach() 
    for (scalar po in pol)
      po[] = 1e-3*noise();
  }

  // if periodic BC, the average of po must vanish
  for (scalar po in pol){
    stats s = statsf(po);
    foreach()
      po[] -= s.sum/s.volume;
  }

  boundary(pol);
  // invert PV at the end of other init event

  if (nptr > 0){
    if ((fp = fopen("ptr0.bas", "r"))) {
      input_matrixl (ptracersl, fp);
      fclose(fp);
    } else {
      foreach() 
        for (scalar ptracers in ptracersl)
          ptracers[] = 1e-3*noise();
    }
    boundary(ptracersl);

    if ((fp = fopen("ptr_relax.bas", "r"))) {
      input_matrixl (ptr_relaxl, fp);
      fclose(fp);
    }
  }

}

/**
   Write parameters
 */
event write_const (t = 0) {
  backup_config();
}

event writestdout (i++) {
/* event writestdout (i=1) { */
  scalar po = pol[0];
  double ke = 0;
  foreach(reduction(+:ke))
    ke -= 0.5*po[]*laplacian(po)*sq(Delta);

  fprintf (stdout,"i = %i, dt = %g, t = %g, ke_1 = %g\n", i, dt, t, ke);
}


event output (t = 0; t <= tend+1e-10;  t += dtout) {
  fprintf(stdout,"write file\n");

  invertq(pol,qol);

  char name[80];
  sprintf (name,"%spo%09d.bas", dpath, i);
  write_field(pol, name, 0.);

  sprintf (name,"%sqo%09d.bas", dpath, i);
  write_field(qol, name, 0.);

  if (dtflt > 0) {
    invertq(tmpl,qofl);
    sprintf (name,"%spf%09d.bas", dpath, i);
    write_field(tmpl, name, 0.);
    nbar = 0; // reset filter average
  }

  /* scalar l[]; */
  /* foreach() */
  /*   l[] = level; */
  /* sprintf (name,"%slevel%09d.dat", dpath, i); */
  /* fp = fopen (name, "w"); */
  /* output_field ({l}, fp); */
  /* fclose(fp); */

  if (ediag>-1){
    double idtout = 1/dtout;
    sprintf (name,"%sde_bf%09d.bas", dpath, i);
    write_field(de_bfl, name, idtout);

    sprintf (name,"%sde_vd%09d.bas", dpath, i);
    write_field(de_vdl, name, idtout);

    sprintf (name,"%sde_j1%09d.bas", dpath, i);
    write_field(de_j1l, name, idtout);

    sprintf (name,"%sde_j2%09d.bas", dpath, i);
    write_field(de_j2l, name, idtout);

    sprintf (name,"%sde_j3%09d.bas", dpath, i);
    write_field(de_j3l, name, idtout);

    sprintf (name,"%sde_ft%09d.bas", dpath, i);
    write_field(de_ftl, name, idtout);

    reset_layer_var(de_bfl);
    reset_layer_var(de_vdl);
    reset_layer_var(de_bfl);
    reset_layer_var(de_j1l);
    reset_layer_var(de_j2l);
    reset_layer_var(de_j3l);
    reset_layer_var(de_ftl);
  }

  if (nptr > 0){
    sprintf (name,"%sptr%09d.bas", dpath, i);
    write_field(ptracersl, name, 0.);
  }

}

/* event adapt (t+=10) { */
/*  astats s = adapt_wavelet (pol, (double []){1e0, 1e0}, maxlevel = 9); */
/*  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc); */
/* } */
