#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re,sys

plt.ion()

dir0 = "../outdir_"

if len(sys.argv) > 1:
  dir0 = dir0 + str(format(sys.argv[1])).zfill(4) + '/'

exec(open(dir0 + "params.in").read())


filep = 'po*'

allfilesp = sorted(glob.glob(dir0 + filep));
nb_files  = len(allfilesp);

p = np.fromfile(allfilesp[0],'f4')
N = int(p[0])
N1 = N + 1
nl = int(len(p)/N1**2)

# grid
Delta = L0/N
x = np.linspace(0.5*Delta, L0-0.5*Delta,N)
y = np.linspace(0.5*Delta, L0-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

# Rossby number (function of latitude)
if (Rom < 0) :
  Ro = 0*yc - Rom
else:
  Ro = Rom/(1 + Rom*beta*(yc-0.5*L0))

# choose iteration (-1: last) and layer (0: upper layer)
nt = 2
il = 0

# load non dimensional pressure p(z,y,x)
p  = np.fromfile(allfilesp[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
# skip first row and first column (coordinates)
p  = p[:,1:,1:]

# horizontal velocity
u = np.gradient(p, axis = 1)/Delta*Ro.reshape(1,N,N)
v = np.gradient(p, axis = 2)/Delta*Ro.reshape(1,N,N)

# kinetic energy
ke = 0.5*(u**2 + v**2)

# buoyancy TODO: missing dh factor
b = np.diff(p,1,0)
# potential energy TODO: missing N^2 factor
pe = b**2

psi = ke[0,:,:]

plt.figure()
plt.imshow(psi[::-1,:])
plt.colorbar()
