/**
QG code

This file is the driver. see qg.h for documentation.

This version uses fields defined on vertices. 
You need to apply the patch "basilisk_bc_vertex.patch" to basilisk.
It is still experimental.

compile with 
qcc -lm -lnetcdf -O3 qg.c -I$DOCUMENT_ROOT/sandbox

 */


#include "grid/multigrid.h"
#include "qg.h"
#include "extra.h"
#include "netcdf_vertex_bas.h"

char* fileout = "vars.nc";

int main(int argc,char* argv[]) {

  // Search for the configuration file with a given path or read params.in 
  if (argc == 2) {
    read_params(argv[1]);
  } else {
    read_params("params.in");
  }

  if (sbc == -1) {
    periodic(right);
    periodic(top);
  }

  create_outdir();

  init_grid (N);
  size(L0);

  run();

/*   set_bc(); */
/*   set_vars(); */
/*   set_const();  */

/*   foreach_vertex()  */
/*     psi[] = 1e-3*noise(); */
/* //    psi[] = 1e-3*noise(); */
/*   boundary({psi}); */
/*   comp_q(psi, q); */
/*   set_bc(); */
/*   boundary({psi}); */
/*   boundary({q}); */
/*   backup_config(); */

/*   sprintf (file_nc,"%s%s", dpath, fileout); */
/*   scalar_list_nc = list_copy({psi, q}); */
/*   create_nc(); */

/*   invertq(psi, q); */

/*   write_nc(); */

}


/**
   Initial conditions
*/
event init (i = 0) {

  foreach_vertex() 
    psi[] = 1e-3*noise();
//    psi[] = 1e-3*noise();
  boundary({psi});

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
  foreach_vertex(reduction(+:ke))
    ke -= 0.5*psi[]*laplacian(psi)*sq(Delta);

  fprintf (stdout,"i = %i, dt = %g, t = %g, ke_1 = %g\n", i, dt, t, ke);
}

event output (t = 0; t <= tend+1e-10;  t += dtout) {
  fprintf(stdout,"write file\n");
  write_nc();
}
