/**
   Netcdf interface for basilisk
*/


//echo | gcc -E -xc -include 'stddef.h' - | grep size_t

#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#pragma autolink -lnetcdf

#if LAYERS
#define NDIMS 4
#else
#define NDIMS 3
#endif

#define Y_NAME "y"
#define X_NAME "x"
#define REC_NAME "time"
#define LVL_NAME "level"


/* /\* For the units attributes. *\/ */
/* #define UNITS "units" */
/* #define PRES_UNITS "hPa" */
/* #define TEMP_UNITS "celsius" */
/* #define MAX_ATT_LEN 80 */

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
int nc_err;
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return;}

// User global variables
scalar * scalar_list_nc;
char file_nc[80];

/* IDs for the netCDF file, dimensions, and variables. */
int ncid;
int t_varid;

// temporary
int nc_varid[1000];
char * nc_varname[1000];
int nvarout = 0;
int nc_rec = -1;

// tempo
int nl_tmp = 1;

void create_nc(scalar * list_out, char* file_out)
{

   // make it global variable
   sprintf (file_nc,"%s", file_out);
   scalar_list_nc = list_copy(list_out);

  if (pid() == 0) { // master

   /* LOCAL IDs for the netCDF file, dimensions, and variables. */
   int x_dimid, y_dimid, lvl_dimid, rec_dimid;
   int y_varid, x_varid, lvl_varid;
   int dimids[NDIMS];

   /* Create the file. */
   if ((nc_err = nc_create(file_nc, NC_CLOBBER, &ncid)))
      ERR(nc_err);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/

   int npx = Dimensions.x;
   int npy = Dimensions.y;
   int Nloc = (1 << depth());
   int Nx = Nloc*npx;
   int Ny = Nloc*npy;


#if LAYERS
   if ((nc_err = nc_def_dim(ncid, LVL_NAME, nl, &lvl_dimid)))
      ERR(nc_err);
#endif
   if ((nc_err = nc_def_dim(ncid, Y_NAME, Ny+1, &y_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, X_NAME, Nx+1, &x_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
      ERR(nc_err);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&y_dimid) and
      similarly for (&x_dimid). */
   if ((nc_err = nc_def_var(ncid, REC_NAME, NC_FLOAT, 1, &rec_dimid,
        		    &t_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, Y_NAME, NC_FLOAT, 1, &y_dimid,
        		    &y_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, X_NAME, NC_FLOAT, 1, &x_dimid,
        		    &x_varid)))
      ERR(nc_err);
#if LAYERS
   if ((nc_err = nc_def_var(ncid, LVL_NAME, NC_FLOAT, 1, &lvl_dimid,
        		    &lvl_varid)))
      ERR(nc_err);
#endif
   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = rec_dimid;
#if LAYERS
   dimids[1] = lvl_dimid;
   dimids[2] = y_dimid;
   dimids[3] = x_dimid;
#else
   dimids[1] = y_dimid;
   dimids[2] = x_dimid;
#endif

   /* Define the netCDF variables */
//   char * str1;
   for (scalar s in scalar_list_nc){
     /* if (strcmp(str1,s.name) != 0) { */
       if ((nc_err = nc_def_var(ncid, s.name, NC_FLOAT, NDIMS,
                                dimids, &nc_varid[nvarout])))
         ERR(nc_err);
     /*   nc_varname[nvarout] = strdup(s.name); */
       nvarout += 1;
     /*   str1 = strdup(s.name); */
     /* } */
   }

   /* /\* Assign units attributes to the netCDF variables. *\/ */
   /* if ((nc_err = nc_put_att_text(ncid, pres_varid, UNITS,  */
   /*      			 strlen(PRES_UNITS), PRES_UNITS))) */
   /*    ERR(nc_err); */
   /* if ((nc_err = nc_put_att_text(ncid, temp_varid, UNITS,  */
   /*      			 strlen(TEMP_UNITS), TEMP_UNITS))) */
   /*    ERR(nc_err); */

   /* End define mode. */
   if ((nc_err = nc_enddef(ncid)))
      ERR(nc_err);

   /*  write coordinates*/
   float yc[Ny+1], xc[Nx+1];
   double Delta = L0*1.0/Nx;
   for (int i = 0; i < Nx+1; i++)
     xc[i] = X0 + i*Delta;

   for (int i = 0; i < Ny+1; i++)
      yc[i] = Y0 + i*Delta;

   if ((nc_err = nc_put_var_float(ncid, y_varid, &yc[0])))
      ERR(nc_err);
   if ((nc_err = nc_put_var_float(ncid, x_varid, &xc[0])))
      ERR(nc_err);

#if LAYERS
   float zc[nl];
   for (int i = 0; i < nl; i++){
      zc[i] = i;
   }
   if ((nc_err = nc_put_var_float(ncid, lvl_varid, &zc[0])))
      ERR(nc_err);
#endif

   /* Close the file. */
   if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
   fprintf(stdout,"*** SUCCESS creating example file %s!\n", file_nc);
  } // end master

}


void write_nc() {

  int npx = Dimensions.x;
  int npy = Dimensions.y;
  int Nloc = (1 << depth());
  int Nx = Nloc*npx;
  int Ny = Nloc*npy;

  int Nx_out = Nx + 1;
  int Ny_out = Ny + 1;

  if (pid() == 0) { // master
    /* open file. */
    if ((nc_err = nc_open(file_nc, NC_WRITE, &ncid)))
      ERR(nc_err);
  }

  // write time
  nc_rec += 1;
  float loctime = t;

  size_t startt[1], countt[1];
  startt[0] = nc_rec; //time
  countt[0] = 1;
  if (pid() == 0) { // master
    if ((nc_err = nc_put_vara_float(ncid, t_varid, startt, countt,
                                    &loctime)))
      ERR(nc_err);
  }



  float * field = (float *)malloc(Nx_out*Ny_out*nl*sizeof(float));

  /* The start and count arrays will tell the netCDF library where to
     write our data. */
  size_t start[NDIMS], count[NDIMS];


  /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
  start[0] = nc_rec; //time
#if LAYERS
  start[1] = 0;     //level
  start[2] = 0;      //y
  start[3] = 0;      //x
#else
  start[1] = 0;      //y
  start[2] = 0;      //x
#endif


  count[0] = 1;
#if LAYERS
  count[1] = nl;
  count[2] = Ny_out;
  count[3] = Nx_out;
#else
  count[1] = Ny_out;
  count[2] = Nx_out;
#endif
  int nv = -1;
  /* char * str1; */
//  foreach_layer() {
  for (scalar s in scalar_list_nc){
    nv += 1;

    for (int k = 0; k < nl; k++) {
      for (int j = 0; j < Ny_out; j++) {
        for (int i = 0; i < Nx_out; i++) {
          field[Ny_out*Nx_out*k + Nx_out*j + i] = nodata;
        }
      }
    }

#if LAYERS
    foreach_layer()
#else
      int _layer = 0;
#endif
      foreach_vertex(cpu){
        int i = (x - X0)/L0*(Nx_out-1);
        int j = (y - Y0)/(L0*(npy*1.0/npx))*(Ny_out-1);
        field[Ny_out*Nx_out*_layer + Nx_out*j + i] = s[];
    }
    if (pid() == 0) { // master
#if _MPI
        MPI_Reduce (MPI_IN_PLACE, &field[0], Ny_out*Nx_out*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif

     if ((nc_err = nc_put_vara_float(ncid, nc_varid[nv], start, count,
        			      &field[0])))
         ERR(nc_err);

  }
#if _MPI
  else // slave
  MPI_Reduce (&field[0], NULL, Ny_out*Nx_out*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif
//  }
  }
//  matrix_free (field);
  free(field);


   /* Close the file. */
  if (pid() == 0) { // master
    if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
  }
//   printf("*** SUCCESS writing example file %s -- %d!\n", file_nc, nc_rec);
}


/**
   Read NC
 */

void read_nc(scalar * list_in, char* file_in, bool read_time){

  int i, ret;
  int ncfile, ndims, nvars, ngatts, unlimited;
  int var_ndims, var_natts;
  int t_id;
  nc_type type;
  char varname[NC_MAX_NAME+1];
  int *dimids = NULL;

  int npx = Dimensions.x;
  int npy = Dimensions.y;
  int Nloc = (1 << depth());
  int Nx = Nloc*npx;
  int Ny = Nloc*npy;



  int Nx_out = Nx+1;
  int Ny_out = Ny+1;

  //  float ** field = matrix_new (Nloc, Nloc, sizeof(float));
  float * field = (float *)malloc(Nx_out*Ny_out*nl*sizeof(float));

  if ((nc_err = nc_open(file_in, NC_NOWRITE, &ncfile)))
    ERR(nc_err);

  if ((nc_err = nc_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited)))
    ERR(nc_err);


  if (read_time) {
    if ((nc_err = nc_inq_varid(ncfile, "time", &t_id)))
      ERR(nc_err);
    size_t startt[1], countt[1];
    startt[0] = 0; //time
    countt[0] = 1;
    float loctime;

    if ((nc_err = nc_get_vara_float(ncfile, t_id, startt, countt,
                                    &loctime)))
      ERR(nc_err);
    t = loctime;
  }

  for (scalar s in list_in){
    for(i=0; i<nvars; i++) {


      if ((nc_err = nc_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                               &var_natts)))
        ERR(nc_err);

      if (strcmp(varname,s.name) == 0) {
        fprintf(stdout,"Reading variable  %s!\n", s.name);

        if (var_ndims < 2) { // error
          fprintf(stdout, "Error: input variables should be at least (x,y)\n");
          return;
        }

        // TODO: Determine if the variable has a time dimension (assume it is unlimited dimension)
        /* int has_time = 0; */
        /* int time_dimid = -1; */
        /* if (nc_inq_unlimdim(ncfile, &time_dimid) == NC_NOERR && time_dimid >= 0) { */
        /*   has_time = (dimids[0] == time_dimid); */
        /* } */

        int has_time = 1;
        // Effective spatial dims = var_ndims - has_time
        int spatial_ndims = var_ndims - has_time;

        if (spatial_ndims < 2) {
          fprintf(stdout, "Error: input variables should have at least (x,y) spatial dims\n");
          return;
        }

        if (spatial_ndims == 2) {
          if (has_time) {
            // (t, y, x)
            size_t start[3] = {0, 0, 0};
            size_t count[3] = {1, Ny_out, Nx_out};
            if ((nc_err = nc_get_vara_float(ncfile, i, start, count, &field[0])))
              ERR(nc_err);
          } else {
            // (y, x)
            size_t start[2] = {0, 0};
            size_t count[2] = {Ny_out, Nx_out};
            if ((nc_err = nc_get_vara_float(ncfile, i, start, count, &field[0])))
              ERR(nc_err);
          }

          foreach_vertex(){
            int i = (x - X0)/L0*(Nx_out - 1);
            int j = (y - Y0)/(L0*(npy*1.0/npx))*(Ny_out - 1);
            s[] = field[Nx_out*j + i];
          }

        } else if (spatial_ndims == 3) {
          if (has_time) {
            // (t, z, y, x)
            size_t start[4] = {0, 0, 0, 0};
            size_t count[4] = {1, nl, Ny_out, Nx_out};
            if ((nc_err = nc_get_vara_float(ncfile, i, start, count, &field[0])))
              ERR(nc_err);
          } else {
            // (z, y, x)
            size_t start[3] = {0, 0, 0};
            size_t count[3] = {nl, Ny_out, Nx_out};
            if ((nc_err = nc_get_vara_float(ncfile, i, start, count, &field[0])))
              ERR(nc_err);
          }

#if LAYERS
          foreach_layer() {
#else
            int _layer = 0;
#endif
            foreach_vertex(){
              int i = (x - X0)/L0*(Nx_out - 1);
              int j = (y - Y0)/(L0*(npy*1.0/npx))*(Ny_out - 1);
              s[] = field[Ny_out*Nx_out*_layer + Nx_out*j + i];
            }
#if LAYERS
          } // foreach_layer
#endif
        } // if spatial_ndims
      } // if (strcmp(varname,s.name) == 0) {
    }
  }

  //  matrix_free (field);
  free(field);

  if ((nc_err = nc_close(ncfile)))
    ERR(nc_err);

  // BD: there should be a call to boundary here (not automatic and needed for MPI)
  boundary(list_in);

}

event cleanup (t = end)
{
  free(scalar_list_nc);
}
