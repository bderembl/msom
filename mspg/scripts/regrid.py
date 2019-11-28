#!/usr/bin/env python

# bi-cubic interpolation of PG fields on fine QG grid

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import sys

# new grid
if len(sys.argv) > 2:
  Nn = int(sys.argv[1])
  nl = int(sys.argv[2])
else:
  print("usage: python regrid.py N nl")
  print("with N the new grid size and nl the number of layer")
  sys.exit(1)

filepo = "psipg_" + str(nl) + "l.bas"
filefr = "frpg_" + str(nl) + "l.bas"
filerd = "rdpg_" + str(nl) + "l.bas"

p = np.fromfile(filepo,'f4')
N = int(p[0])
N1 = N + 1
N2 = N + 2
nl = int(len(p)/N1**2)

# old grid
Delta = 1/N
x  = np.linspace(0.5*Delta, 1-0.5*Delta,N)
x2 = np.linspace(-0.5*Delta, 1+0.5*Delta,N2)

# new grid
Deltan = 1/Nn
xn = np.linspace(0.5*Deltan, 1-0.5*Deltan,Nn)

po = np.fromfile(filepo,'f4').reshape(nl,N1,N1).transpose(0,2,1)
fr = np.fromfile(filefr,'f4').reshape(nl,N1,N1).transpose(0,2,1)
rd = np.fromfile(filerd,'f4').reshape(N1,N1).transpose(1,0)

po2 = np.zeros((nl,N2,N2))
fr2 = np.zeros((nl,N2,N2))
rd2 = np.zeros((N2,N2))

po2[:,1:-1,1:-1] = po[:,1:,1:]
fr2[:,1:-1,1:-1] = fr[:,1:,1:]
rd2[1:-1,1:-1] = rd[1:,1:]

# boundaries

po2[:,0,:]  = -po2[:,1,:]
po2[:,-1,:] = -po2[:,-2,:]
po2[:,:,0]  = -po2[:,:,1]
po2[:,:,-1] = -po2[:,:,-2]

fr2[:,0,:]  = fr2[:,1,:]
fr2[:,-1,:] = fr2[:,-2,:]
fr2[:,:,0]  = fr2[:,:,1]
fr2[:,:,-1] = fr2[:,:,-2]

rd2[0,:]  = rd2[1,:]
rd2[-1,:] = rd2[-2,:]
rd2[:,0]  = rd2[:,1]
rd2[:,-1] = rd2[:,-2]

# corners
po2[:,0,0]   = -po2[:,0,1]   - po2[:,1,0]   - po2[:,1,1]
po2[:,-1,0]  = -po2[:,-1,1]  - po2[:,-2,0]  - po2[:,-2,1]
po2[:,0,-1]  = -po2[:,1,-1]  - po2[:,0,-2]  - po2[:,1,-2]
po2[:,-1,-1] = -po2[:,-1,-2] - po2[:,-2,-2] - po2[:,-2,-1]

fr2[:,0,0]   = fr2[:,1,1]
fr2[:,-1,0]  = fr2[:,-2,1]
fr2[:,0,-1]  = fr2[:,1,-2]
fr2[:,-1,-1] = fr2[:,-2,-2]

rd2[0,0]   = rd2[1,1]
rd2[-1,0]  = rd2[-2,1]
rd2[0,-1]  = rd2[1,-2]
rd2[-1,-1] = rd2[-2,-2]

# interpolate
po3 = np.zeros((nl,Nn,Nn))
fr3 = np.zeros((nl,Nn,Nn))
rd3 = np.zeros((Nn,Nn))

for il in range(0,nl):
  interp_spline = RectBivariateSpline(x2, x2, po2[il,:,:])
  po3[il,:,:] = interp_spline(xn, xn)

  interp_spline = RectBivariateSpline(x2, x2, fr2[il,:,:])
  fr3[il,:,:] = interp_spline(xn, xn)

interp_spline = RectBivariateSpline(x2, x2, rd2[:,:])
rd3[:,:] = interp_spline(xn, xn)

#save output
Fr_o = np.zeros((nl,Nn+1,Nn+1))
Fr_o[:,1:,1:] = fr3
Fr_o[:,0,:] = 0
Fr_o[:,:,0] = 0
Fr_o[:,0,0] = Nn
Fr_o = np.transpose(Fr_o,(0,2,1))
fileFr = 'frpg_' + str(nl) +'l_N' + str(Nn) + '.bas'
Fr_o.astype('f4').tofile(fileFr)

Rd_o = np.zeros((Nn+1,Nn+1))
Rd_o[1:,1:] = rd3
Rd_o[0,:] = 0
Rd_o[:,0] = 0
Rd_o[0,0] = Nn
Rd_o = np.transpose(Rd_o,(1,0))
fileRd = 'rdpg_' + str(nl) +'l_N' + str(Nn) + '.bas'
Rd_o.astype('f4').tofile(fileRd)

psi_o = np.zeros((nl,Nn+1,Nn+1))
psi_o[:,1:,1:] = po3
psi_o[:,0,:] = 0
psi_o[:,:,0] = 0
psi_o[:,0,0] = Nn
psi_o = np.transpose(psi_o,(0,2,1))
fileppg = 'psipg_' + str(nl) +'l_N' + str(Nn) + '.bas'
psi_o.astype('f4').tofile(fileppg)
