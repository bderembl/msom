#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import qg as bas

dir0 = "./"

params_file = "params.in"

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

# Adams Bashforth arrays
si_t = int(tend/DT)
F1 = np.zeros((nl,N,N))
F2 = np.zeros((nl,N,N))
F3 = np.zeros((nl,N,N))

# forward integration: direction = 1.0
# backward integration: direction = -1.0
direction = 1.0

# time integration loop
for nt in range(0,si_t):
  
  # compute RHS
  bas.pystep_bfn(p, F1, direction)

  # BNF nudging goes here
  # F1 = F1 +  nudging

  # Adams Bashforth time stepping
  p = p + DT/12*(23*F1 - 16*F2 + 5*F3)
  F3 = np.copy(F2)
  F2 = np.copy(F1)


# the following two lines will cause a segmentation fault if the
# script is run interactively

bas.trash_vars()
bas.trash_vars_bfn()
