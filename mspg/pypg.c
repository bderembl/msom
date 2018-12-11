/**
   Planetary geostrophic model

   Python driver
*/

#include "grid/multigrid.h"
#include "../msqg/auxiliar_input.h"
#include "pg.h"


double k (double x, double y, double s) { return (5e-4);}

double tau0 = 0.12;

double taux   (double x, double y){ return (tau0*sin(2*(y-ys)*pi));}
double taux_y (double x, double y){ return (2*pi*tau0*cos(2*(y-ys)*pi));}
double tauy   (double x, double y){ return (0.);}
double tauy_x (double x, double y){ return (0.);}

