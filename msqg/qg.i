%module qg
%include "common.i"
%include "grid/multigrid.i"
%include "predictor-corrector.i"
%include "poisson.i"
%include "timestep.i"


%apply (double * IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * val1, int len1, int len2, int len3)
}

%apply (double * INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  (double * val2, int len4, int len5, int len6)
}

%inline %{
  void pystep ( double * val1, int len1, int len2, int len3,
                double * val2, int len4, int len5, int len6);
%}
%{
  extern scalar sig_lev;
  extern scalar sig_filt;
  extern scalar Rd;
  extern scalar Ro;
%}

extern scalar sig_lev;
extern scalar sig_filt;
extern scalar Rd;
extern scalar Ro;

%pythoncode %{
sig_lev = scalar(_qg.cvar.sig_lev)
sig_filt = scalar(_qg.cvar.sig_filt)
Rd = scalar(_qg.cvar.Rd)
Ro = scalar(_qg.cvar.Ro)
%}

%inline %{
  void set_vars ();
  void set_const ();
  void read_params();
  void create_outdir();
  void backup_config();
  void trash_vars();
%}

