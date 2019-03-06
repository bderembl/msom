
/**
   Interpolation routine

   qcc -lm -O3 regrid.c -o regrid.e
   ./regrid.e 512 0 psipg.bas
   ./regrid.e 512 1 frpg.bas
*/

#include "run.h"
#include "../../msqg/auxiliar_input.h"
#include <stdio.h>

scalar * omegal;
char name[80];
int Nout; 
int ibc = 1; 

int main(int argc, char **argv) {
  
  if (argc >= 4) {
    Nout = atoi(argv[1]); 
    ibc = atoi(argv[2]); 
    strcpy(name, argv[3]);
    } else {
      fprintf(stderr, "Usage: ./regrid.e N_output bc filename.\n");
      fprintf(stderr, "with bc = 0: dirichlet BC, 1, neumann\n");
      exit(0);
    }

  FILE * fp = fopen (name, "r"); 
  float nsize;
  fread(&nsize, sizeof(float), 1, fp);
  fseek(fp, 0, SEEK_END); // seek to end of file
  int nl2 = ftell(fp); // get current file pointer
  fseek(fp, 0, SEEK_SET); // seek back to beginning of file
  fclose(fp);

  N = (int) nsize;
  int nl = nl2/sizeof(float)/N/N;
  printf("Input field of size N = %d, nl= %d\n", N, nl);


  init_grid (N);

  for (int l = 0; l < nl; l++) {
    scalar omega = new scalar;
    omegal = list_append (omegal, omega);
    if (ibc == 0) {
      omega[top]    = dirichlet(0);
      omega[bottom] = dirichlet(0);
      omega[right]  = dirichlet(0);
      omega[left]   = dirichlet(0);
    }
  }

  fp = fopen (name, "r");  
  input_matrixl (omegal, fp);
  fclose(fp);

  boundary(omegal);

  char name2[80];
  int len = strlen(name);
  name[len - 4] = '\0'; // strip .bas
  sprintf (name2,"%s_N%d.bas", name, Nout);

  fp = fopen (name2, "w");
  output_matrixl (omegal, fp, Nout, linear=1);
  fclose(fp);

  
}
