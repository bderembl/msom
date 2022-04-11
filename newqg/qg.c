/**
QG code

compile with 

qcc -lm -lnetcdf -O3 qg.c

 */


#include "grid/multigrid.h"
#include "qg.h"
#include "extra.h"
#include "netcdf_bas.h"
#include "auxiliar_io.h"

char* fileout = "vars.nc";

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
    input_matrixl (psi, fp);
    fclose(fp);
  } else {
  foreach() 
    foreach_layer() 
      psi[] = 1e-3*noise();
  }

}

/**
   Write parameters
 */
event write_const (t = 0) {
  backup_config();

  sprintf (file_nc,"%s%s", dpath, fileout);
  scalar_list_nc = list_copy({psi, q});
  create_nc();
}

event writestdout (i++) {
  double ke = 0;
  foreach(reduction(+:ke))
    ke -= 0.5*psi[]*laplacian(psi)*sq(Delta);

  fprintf (stdout,"i = %i, dt = %g, t = %g, ke_1 = %g\n", i, dt, t, ke);
}

event output (t = 0; t <= tend+1e-10;  t += dtout) {
  fprintf(stdout,"write file\n");
  write_nc();
}
