/**
   ## Layerered variables initialization, etc
*/

scalar * create_layer_var (scalar * psil, int nl, int bc_type)
{
  assert (psil == NULL);
  assert (nl > 0);

  for (int l = 0; l < nl; l++) {
    scalar po = new scalar;
    psil = list_append (psil, po);
    /**
       By convention, psi = 0 at the boundary which is located 1/2 grid
       point outside the computational domain
     */
    if (bc_type == 0){
      po[right]  = dirichlet(0);
      po[left]   = dirichlet(0);
      po[top]    = dirichlet(0);
      po[bottom] = dirichlet(0);
      // instead of po[...] = 0.; -> boundary located 1 grid point outside
      /* po[right]  = 0; */
      /* po[left]   = 0; */
      /* po[top]    = 0; */
      /* po[bottom] = 0; */
    }
  }

  foreach() 
    for (scalar po in psil) {po[] = 0.0;} 
  boundary(psil);

  return psil;
}

void reset_layer_var(scalar *psil)
{
  foreach()
    for (scalar po in psil) {po[] = 0.0;}
}


void list_copy_deep (scalar * listin, scalar * listout, int nl)
{
  foreach()
    for (int l = 0; l < nl ; l++) {
      scalar listi = listin[l];
      scalar listo = listout[l];
      listo[] = listi[];
    }
}
