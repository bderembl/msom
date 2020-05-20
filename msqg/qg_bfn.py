#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import qg as bas

plt.ion()

dir0 = "./"

params_file = "params.in"

flag_q = 1  # 0: integration in p, 1: integration in q

# read parameters in python
exec(open(dir0 + params_file).read())

# init basilisk variables
bas.read_params(params_file)
bas.init_grid(N)
bas.set_vars()
bas.set_vars_bfn()
bas.set_const()

# local grid
Delta = L0/N
x = np.linspace(0.5*Delta,L0 - 0.5*Delta,N)
xc,yc = np.meshgrid(x,x)

# initial condition
p  = np.zeros((nl,N,N))
p[0,:,:] = np.sin(2*np.pi/L0*xc)*np.sin(2*np.pi/L0*yc)

q  = np.zeros((nl,N,N))
if flag_q:
  bas.pyp2q(p,q)

# Adams Bashforth arrays
si_t = int(tend/DT)
F1 = np.zeros((nl,N,N))
F2 = np.zeros((nl,N,N))
F3 = np.zeros((nl,N,N))
if flag_q:
  var = np.copy(q)
else:
  var = np.copy(p)

# forward integration: direction = 1.0
# backward integration: direction = -1.0
direction = 1.0

# time integration loop
for nt in range(0,10):
  
  # compute RHS
  bas.pystep_bfn(var, F1, direction, flag_q)

  # BNF nudging goes here
  # F1 = F1 +  nudging

  # Adams Bashforth time stepping
  var = var + DT/12*(23*F1 - 16*F2 + 5*F3)
  F3 = np.copy(F2)
  F2 = np.copy(F1)

if flag_q:
  q = np.copy(var)
  bas.pyq2p(p,q)
else:
  p = np.copy(var)
  bas.pyp2q(p,q)

# the following two lines will cause a segmentation fault if the
# script is run interactively

bas.trash_vars()
bas.trash_vars_bfn()
