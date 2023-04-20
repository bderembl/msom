/**
QG code

This file is the driver. see qg.h for documentation.

This version uses fields defined on vertices. 
It is still experimental.

compile with 
qcc -lm -lnetcdf -O3 qg.c -I$DOCUMENT_ROOT/sandbox

create a restart file:
ncks -d time,198,198 vars.nc restart.nc


 */

#if LAYERS
// condition
#else
int nl = 1;
#endif



#include "grid/multigrid.h"
#include "extra.h"
#include "netcdf_vertex_bas.h"

#include "qg.h"
//#include "qg_barotropic.h"
#include "qg_baroclinic_ms.h"

char* fileout = "vars.nc";

int main(int argc,char* argv[]) {


  // declare user parameters
  params = array_new();
  add_param ("N", &N, "int");
  add_param ("nl", &nl, "int");
  add_param ("flag_ms", &flag_ms, "int");
  add_param ("L0", &L0, "double");
  add_param ("f0", &f0, "double");
  add_param ("beta", &beta, "double");
  add_param ("nu", &nu, "double");
  add_param ("hEkb", &hEkb, "double");
  add_param ("tau0", &tau0, "double");
  add_param ("noise_init", &noise_init, "double");
  add_param ("Lfmax", &Lfmax, "double");
  add_param ("dtflt", &dtflt, "double");
  add_param ("dh", &dh[0], "array");
  add_param ("N2", &N2[0], "array");
  add_param ("bc_fac", &bc_fac, "double");
  add_param ("DT", &DT, "double");
  add_param ("tend", &tend, "double");
  add_param ("dtout", &dtout, "double");
  add_param ("CFL", &CFL, "double");
  add_param ("TOLERANCE", &TOLERANCE, "double");
  add_param ("dtdiag", &dtdiag, "double");


  // Search for the configuration file with a given path or read params.in 
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in

  read_params(file_param);
  create_outdir();
  backup_config(file_param);


  if (bc_fac == -1) {
    periodic(right);
    periodic(top);
  }
 
  init_grid (N);
  size(L0);

  run();

  array_free (params);

}


/**
   Surface forcing
*/
event init (i = 0) {

  foreach_vertex()
    q_forcing[] = -tau0/dh[0]*2*pi/L0*sin(2*pi*y/L0);
//    q_forcing[] = -tau0/L0*pi*sin(pi*y/L0);

}

/**
   Write parameters
 */
event write_const (t = 0) {

  sprintf (file_nc,"%s%s", dpath, fileout);
  scalar_list_nc = list_copy({psi, q});
  create_nc();
}

event output (t = 0; t <= tend+1e-10;  t += dtout) {
  fprintf(stdout,"write file\n");

  if (i == 0)
    invert_q(psi, q);

  write_nc();

  fprintf(stdout,"file written \n");
}

event writestdout (i++) {
  double ke = 0;
  foreach_vertex(reduction(+:ke))
    ke -= 0.5*psi[]*laplacian(psi)*sq(Delta);

  fprintf (stdout,"i = %i, dt = %g, t = %g, ke_1 = %g\n", i, dt, t, ke);
}

