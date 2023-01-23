/**
QG code

This file is the driver. see qg.h for documentation.

This version uses fields defined on vertices. 
You need to apply the patch "basilisk_bc_vertex.patch" to basilisk.
It is still experimental.

compile with 
qcc -lm -lnetcdf -O3 qg.c -I$DOCUMENT_ROOT/sandbox

create a restart file:
ncks -d time,198,198 vars.nc restart.nc


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
}


/**
   Initial conditions, surface forcing and PG fields
*/
event init (i = 0) {

  foreach_vertex() 
    foreach_layer()
      psi[] = 1e-3*noise();
//      psi[] = 1e-3*noise();

  foreach_vertex()
    q_forcing[] = tau0/dh[0]*3/2*pi/L0*sin(2*pi*y/L0)*sin(pi*y/L0);


  fprintf(stdout, "Read input files:\n");

  FILE * fp;
  char name[80];
  sprintf (name,"restart.nc");
  if ((fp = fopen(name, "r"))) {
    read_nc({psi}, name);
    fclose(fp);
    backup_file(name);
    fprintf(stdout, "%s .. ok\n", name);
  }

  sprintf (name,"psipg_%dl_N%d.nc", nl,N);
  if ((fp = fopen(name, "r"))) {
    read_nc({psi_pg}, name);
    fclose(fp);
    backup_file(name);
    fprintf(stdout, "%s .. ok\n", name);
  }

  sprintf (name,"gp_l_N%d.nc",N);
  if ((fp = fopen(name, "r"))) {
    read_nc({gp_l}, name);
    fclose(fp);
    backup_file(name);
    fprintf(stdout, "%s .. ok\n", name);
  }

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
  invertq(psi, q);
  write_nc();
}
