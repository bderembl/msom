/**
QG code

This file is the driver. see qg.h for documentation.

This version uses fields defined on vertices. 
It is still experimental.

for now compile with -disable-dimensions

compile with 
qcc -lm -lnetcdf -O3 qg.c -o qg.e (-DLAYERS=1) (-fopenmp) (-I$DOCUMENT_ROOT/sandbox) (-disable-dimensions )
export OMP_NUM_THREADS=20 (?)


MPI:
CC99='mpicc -std=c99' qcc -D_MPI=1 -lm -lnetcdf -O3 qg.c -o qg.e -grid=multigrid (-DLAYERS=1)
mpirun -np 16 ./qg.e


HPC:
qcc -D_MPI=1 -grid=multigrid -source qg.c    ( -DLAYERS=1 )
rsync _qg.c
mpicc -Wall -std=c99 -O2 _qg.c -lm -lnetcdf -o qg.e 


create a restart file:
ncks -d time,-1,-1 vars.nc restart.nc


exerimental flags:
-DFORCING_3D=1 to add a 3d forcing (field read in the input_vars...nc file)

 */

#if LAYERS
// condition
#else
int nl = 1;
#endif



#include "grid/multigrid.h"
#include "extra.h"
#include "netcdf_vertex_bas.h"
//#include "pnetcdf_vertex_bas.h"

#include "qg_stochastic.h"
#include "qg.h"

#if LAYERS
#include "qg_baroclinic_ms.h"
#else
#include "qg_barotropic.h"
#endif

char* fileout = "vars.nc";

double tau1 = 0.; //forcing parameter
double tf1 = 1.; //forcing parameter
double tf2 = 1.; //forcing parameter
double dy_ws = 1.; //forcing parameter
double forc_mode = 2.0; //forcing parameter



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
  add_param ("nu4", &nu4, "double");
  add_param ("hEkb", &hEkb, "double");
  add_param ("gp_low", &gp_low, "double");
  add_param ("scale_topo", &scale_topo, "double");
  add_param ("tau0", &tau0, "double");
  add_param ("tau1", &tau1, "double");
  add_param ("tf1", &tf1, "double");
  add_param ("tf2", &tf2, "double");
  add_param ("dy_ws", &dy_ws, "double");
  add_param ("forc_mode", &forc_mode, "double");
  add_param ("noise_init", &noise_init, "double");
  add_param ("Lfmax", &Lfmax, "double");
  add_param ("Lfmin", &Lfmin, "double");
  add_param ("fac_filt_Rd", &fac_filt_Rd, "double");
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
#ifdef _STOCHASTIC
  add_param ("amp_stoch", &amp_stoch, "double");
  add_param ("L_filt", &L_filt, "double");
#endif

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
event forcing (i++) {
//event init (i = 0) {

  foreach_vertex()
//    q_forcing[] = -tau0/L0*pi*sin(pi*y/L0);
//    q_forcing[] = - (tau0 + tau1*cos(2*pi*t/tf1))/dh[0]*2*pi/L0*sin(2*pi*y/L0);
    q_forcing[] = -(tau0 + tau1*cos(2*pi*t/tf1))/dh[0]*forc_mode*pi/L0         \
    *sin(forc_mode*pi*(y + y*(y-L0)*2/(L0*L0)*dy_ws*sin(2*pi*t/tf2))/L0);

}

/**
   Write parameters
 */
event write_const (t = 0) {

  char file_tmp[90];
  sprintf (file_tmp,"%s%s", dpath, fileout);
//  scalar_list_nc = list_copy();
//  scalar_list_nc = list_copy({psi, q, psi_f});
  create_nc({psi, q}, file_tmp);
}

event output (t = 0; t <= tend+1e-10;  t += dtout) {
  fprintf(stdout,"write file\n");

  if (i == 0){
//    wavelet_filter ( q, psi);
    invert_q(psi, q);
  }

  write_nc();
  nbar = 0;

  fprintf(stdout,"file written \n");
}

event writestdout (i++) {
  double ke = 0;
  foreach_vertex(reduction(+:ke))
    ke -= 0.5*psi[]*laplacian(psi)*sq(Delta);

  fprintf (stdout,"i = %i, dt = %g, t = %g, ke_1 = %g\n", i, dt, t, ke);
}

