/**
# Multiple scale quasi geostrophic model

This is the driver file for qg.h.  We define the grid, the topography,
the initial conditions and the output routines.

compile with (openmp)
qcc -lm -llapacke qg.c -O3 -o qg.e -fopenmp
export OMP_NUM_THREADS=20 (?)
./qg.e

MPI:
CC99='mpicc -std=c99' qcc -D_MPI=1 -lm -O3 -llapacke qg.c -o qg.e -grid=multigrid
mpirun -np 16 ./qg.e

HPC:
qcc -D_MPI=1 -D MKL -grid=multigrid -source qg.c
mpicc -Wall -std=c99 -O2 _qg.c -lm -mkl -o qg.e
*/

#include "grid/multigrid.h"
#include "auxiliar_input.h"
#include "qg.h"

// for mkdir
#include <sys/stat.h>
#include <sys/types.h>

int main() {

/**
   Read input parameters
 */

  FILE * fp;
  if (fp = fopen("params.in", "rt")) {
    char tempbuff[100];
    char tmps1[16];
    char tmps2[16];

    while(fgets(tempbuff,100,fp)) {
      sscanf(tempbuff, "%15s = %15s", tmps1, tmps2);
      if      (strcmp(tmps1,"N")    ==0) { N     = atoi(tmps2); }
      else if (strcmp(tmps1,"nl")   ==0) { nl    = atoi(tmps2); }
      else if (strcmp(tmps1,"ediag")==0) { ediag = atoi(tmps2); }
      else if (strcmp(tmps1,"L0")   ==0) { L0    = atof(tmps2); }
      else if (strcmp(tmps1,"Rom")  ==0) { Rom   = atof(tmps2); }
      else if (strcmp(tmps1,"Ek")   ==0) { Ek    = atof(tmps2); }
      else if (strcmp(tmps1,"Re")   ==0) { Re    = atof(tmps2); }
      else if (strcmp(tmps1,"Re4")  ==0) { Re4   = atof(tmps2); }
      else if (strcmp(tmps1,"beta") ==0) { beta  = atof(tmps2); }
      else if (strcmp(tmps1,"afilt")==0) { afilt = atof(tmps2); }
      else if (strcmp(tmps1,"Lfmax")==0) { Lfmax = atof(tmps2); }
      else if (strcmp(tmps1,"DT")   ==0) { DT    = atof(tmps2); }
      else if (strcmp(tmps1,"tend") ==0) { tend  = atof(tmps2); }
      else if (strcmp(tmps1,"dtout")==0) { dtout = atof(tmps2); }
      else if (strcmp(tmps1,"dtflt")==0) { dtflt = atof(tmps2); }
      else if (strcmp(tmps1,"CFL")  ==0) { CFL   = atof(tmps2); }
    }
    fclose(fp);
  } else {
    fprintf(stdout, "file params.in not found\n");
    exit(0);
  }

  /**
     Create output directory and copy input parameter file for backup
  */

  char ch;
  char name[80];
@if _MPI 
  sprintf(dpath, "outdir_%04d/", 1);
  int res = mkdir(dpath, 0777);
@else
  for (int i=1; i<10000; i++) {
    sprintf(dpath, "outdir_%04d/", i);
    if (mkdir(dpath, 0777) == 0) {
      fprintf(stdout,"Writing output in %s\n",dpath);
      break;
    }
  }
@endif

  sprintf (name,"%sparams.in", dpath);
  FILE * source = fopen("params.in", "r");
  FILE * target = fopen(name, "w");
  while ((ch = fgetc(source)) != EOF)
    fputc(ch, target);
  fclose(source);
  fclose(target);

  fprintf(stdout, "Config: N = %d, nl = %d, L0 = %g\n", N, nl, L0);

  init_grid (N);
  size(L0);
  run();
}


event init (i = 0) {
/**
   Layer thickness and large scale variables
*/
  fprintf(stdout, "Read input files:\n");

  char name[80];
  sprintf (name,"dh_%dl.bin", nl);
  float dh[nl];
  FILE * fp = fopen (name, "r");
  fread(&dh, sizeof(float), nl, fp);
  fclose(fp);

  for (int l = 0; l < nl ; l++)
    dhf[l] = dh[l];

  fprintf(stdout, "%s .. ok\n", name);

  sprintf (name,"psipg_%dl_N%d.bas", nl,N);
  fp = fopen (name, "r");
  input_matrixl (ppl, fp);
  fclose(fp);

  fprintf(stdout, "%s .. ok\n", name);

  sprintf (name,"frpg_%dl_N%d.bas", nl,N);
  fp = fopen (name, "r");
  input_matrixl (Frl, fp);
  fclose(fp);

  fprintf(stdout, "%s .. ok\n", name);

  // copy input fields for backup
  sprintf (name,"%sdh_%dl.bin", dpath, nl);
  fp = fopen (name, "w");
  fwrite(&dh, sizeof(float), nl, fp);
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

event filter (t = 0; t <= tend+1e-10;  t += dtflt) {
  fprintf(stdout,"Filter solution\n");
  if (ediag>-1) 
    filter_de (qol, pol, de_ftl);
  wavelet_filter ( qol, pol, qofl, dtflt, nbar)
}

event writestdout (i++) {
/* event writestdout (i=1) { */
  scalar po = pol[0];
  scalar zeta = zetal[0];
  double ke = 0;
  foreach(reduction(+:ke))
    ke -= 0.5*po[]*zeta[]*sq(Delta);

  fprintf (stdout,"i = %i, dt = %g, t = %g, ke_1 = %g\n", i, dt, t, ke);
}

/**
   Write parameters
 */
event write_const (t = 0) {
  char name[80];
  sprintf (name,"%siBu.bas", dpath);
  FILE * fp = fopen (name, "w");
  output_matrixl (iBul, fp);
  fclose(fp);

  sprintf (name,"%ssig_filt.bas", dpath);
  fp = fopen (name, "w");
  output_matrixl ({sig_filt}, fp);
  fclose(fp);

  // copy input field for backup
  sprintf (name,"%spsipg_%dl_N%d.bas", dpath, nl,N);
  fp = fopen (name, "w");
  output_matrixl (ppl, fp);
  fclose(fp);

  sprintf (name,"%sfrpg_%dl_N%d.bas", dpath, nl,N);
  fp = fopen (name, "w");
  output_matrixl (Frl, fp);
  fclose(fp);

}

event output (t = 0; t <= tend+1e-10;  t += dtout) {
  fprintf(stdout,"write file\n");

  invertq(pol,qol);

  char name[80];
  sprintf (name,"%spo%09d.bas", dpath, i);
  write_field(pol, name, 0.);

  sprintf (name,"%sqo%09d.bas", dpath, i);
  write_field(qol, name, 0.);

  invertq(tmpl,qofl);
  sprintf (name,"%spf%09d.bas", dpath, i);
  write_field(tmpl, name, 0.);
  
  nbar = 0; // reset filter average

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
}

/* event adapt (t+=10) { */
/*  astats s = adapt_wavelet (pol, (double []){1e0, 1e0}, maxlevel = 9); */
/*  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc); */
/* } */
