#!/usr/bin/env python

# bi-cubic interpolation of PG fields on fine QG grid

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

# new grid
Nn = 512

filepo = "psipg_4l.bas"
filefr = "frpg_4l.bas"

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

po2 = np.zeros((nl,N2,N2))
fr2 = np.zeros((nl,N2,N2))

po2[:,1:-1,1:-1] = po[:,1:,1:]
fr2[:,1:-1,1:-1] = fr[:,1:,1:]

# boundaries

po2[:,0,:]  = -po2[:,1,:]
po2[:,-1,:] = -po2[:,-2,:]
po2[:,:,0]  = -po2[:,:,1]
po2[:,:,-1] = -po2[:,:,-2]

fr2[:,0,:]  = fr2[:,1,:]
fr2[:,-1,:] = fr2[:,-2,:]
fr2[:,:,0]  = fr2[:,:,1]
fr2[:,:,-1] = fr2[:,:,-2]

# corners
po2[:,0,0]   = -po2[:,0,1]   - po2[:,1,0]   - po2[:,1,1]
po2[:,-1,0]  = -po2[:,-1,1]  - po2[:,-2,0]  - po2[:,-2,1]
po2[:,0,-1]  = -po2[:,1,-1]  - po2[:,0,-2]  - po2[:,1,-2]
po2[:,-1,-1] = -po2[:,-1,-2] - po2[:,-2,-2] - po2[:,-2,-1]

fr2[:,0,0]   = fr2[:,1,1]
fr2[:,-1,0]  = fr2[:,-2,1]
fr2[:,0,-1]  = fr2[:,1,-2]
fr2[:,-1,-1] = fr2[:,-2,-2]

# interpolate
po3 = np.zeros((nl,Nn,Nn))
fr3 = np.zeros((nl,Nn,Nn))

for il in range(0,nl):
  interp_spline = RectBivariateSpline(x2, x2, po2[il,:,:])
  po3[il,:,:] = interp_spline(xn, xn)

  interp_spline = RectBivariateSpline(x2, x2, fr2[il,:,:])
  fr3[il,:,:] = interp_spline(xn, xn)

#save output
Fr_o = np.zeros((nl,Nn+1,Nn+1))
Fr_o[:,1:,1:] = fr3
Fr_o[:,0,:] = 0
Fr_o[:,:,0] = 0
Fr_o[:,0,0] = Nn
Fr_o = np.transpose(Fr_o,(0,2,1))
fileFr = 'frpg_' + str(nl) +'l_N' + str(Nn) + '.bas'
Fr_o.astype('f4').tofile(fileFr)

psi_o = np.zeros((nl,Nn+1,Nn+1))
psi_o[:,1:,1:] = po3
psi_o[:,0,:] = 0
psi_o[:,:,0] = 0
psi_o[:,0,0] = Nn
psi_o = np.transpose(psi_o,(0,2,1))
fileppg = 'psipg_' + str(nl) +'l_N' + str(Nn) + '.bas'
psi_o.astype('f4').tofile(fileppg)
