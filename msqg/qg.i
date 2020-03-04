%module qg
%include "common.i"
%include "grid/multigrid.i"
%include "predictor-corrector.i"
%include "poisson.i"
%include "timestep.i"
%include "qg_energy.i"
%include "qg_bfn.i"

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
  void read_params(char* path2file);
  void create_outdir();
  void backup_config();
  void trash_vars();
%}

