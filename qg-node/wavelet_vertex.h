/**
   TODO wavelet version
 */



scalar mask_c[];


void wavelet_mask (scalar s, scalar w)
{
  restriction ({s});
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    foreach_coarse_level (l) {
      foreach_child()
        w[] = s[];
      s.prolongation (point, s);
      foreach_child() {
        double sp = s[];
        s[] = w[];
        /* difference between fine value and its prolongation */
        w[] = (w[] - sp)*mask_c[];
      }
    }
    boundary_level ({w}, l + 1);
  }
  /* root cell */
  foreach_level(0) 
    w[] = s[]*mask_c[];
  boundary_level ({w}, 0);
}

void inverse_wavelet_mask (scalar s, scalar w)
{
  foreach_level(0) 
    s[] = w[]*mask_c[];
  boundary_level ({s}, 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    foreach_coarse_level (l) {
      s.prolongation (point, s);
      foreach_child()
        s[] = (s[] + w[])*mask_c[];
    }
    boundary_level ({s}, l + 1);
  }
}


/* void wavelet_vertex (scalar s, scalar w) */
/* { */
/*   restriction ({s}); */
/*   for (int l = grid->maxdepth - 1; l >= 0; l--) { */
/*     foreach_coarse_vertex_level (l) { */
/*       foreach_child() */
/*         w[] = s[]; */
/*       s.prolongation (point, s); */
/*       foreach_child() { */
/*         double sp = s[]; */
/*         s[] = w[]; */
/*         /\* difference between fine value and its prolongation *\/ */
/*         w[] -= sp; */
/*       } */
/*     } */
/*     boundary_level ({w}, l + 1); */
/*   } */
/*   /\* root cell *\/ */
/*   foreach_vertex_level(0) */
/*     w[] = s[]; */
/*   boundary_level ({w}, 0); */
/* } */

/* void inverse_wavelet_vertex (scalar s, scalar w) */
/* { */
/*   foreach_vertex_level(0) */
/*     s[] = w[]; */
/*   boundary_level ({s}, 0); */
/*   for (int l = 0; l <= grid->maxdepth - 1; l++) { */
/*     foreach_coarse_vertex_level (l) { */
/*       s.prolongation (point, s); */
/*       foreach_child() */
/*         s[] += w[]; */
/*     } */
/*     boundary_level ({s}, l + 1); */
/*   } */
/* } */


/* trace */
/* void wavelet_filter_vertex(scalar q, scalar psi) */
/* { */

/*   invert_q(psi,q); */

/*   foreach_layer(){ */
/*     vertex scalar w[]; */

/*     w[top]    = 0; */
/*     w[bottom] = 0; */
/*     w[right]  = 0; */
/*     w[left]   = 0; */

/* //    psi.restriction = restriction_coarsen_vert2; */

/*     wavelet_vertex(psi,w); */

/*   for (int l = 0; l < 10; l++) { */
/*     foreach_vertex_level (l) */
/*       w[] = 0.; */
/*     boundary_level ({w}, l); */
/*   } */


/*     for (int l = 0; l <= depth(); l++) { */
/*       foreach_vertex_level (l) */
/*         w[] *= sig_lev[]; */
/*       boundary_level ({w}, l); */
/*     } */
/*     inverse_wavelet_vertex (psi, w); */

/* //    psi.restriction = restriction_vert; */
/*   } */
  
/*   comp_q(psi,q); */
  
/* } */
