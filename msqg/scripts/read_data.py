#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re

plt.ion()

dir0 = "./"
filep = 'po*'

# Non dimensional parameters
Lt = 100  # size of the domain (adim)
Rom = 0.025  # rossby number in the middle of the domain
beta = 0.5   # beta (adim)

# reference length scale and velocity scale
lref = 50e3 # m
uref = 0.1  # m/s

allfilesp = sorted(glob.glob(dir0 + filep));
nb_files  = len(allfilesp);

p = np.fromfile(allfilesp[0],'f4')
N = int(p[0])
N1 = N + 1
nl = int(len(p)/N1**2)

# grid
Delta = Lt/N
x = np.linspace(0.5*Delta, Lt-0.5*Delta,N)
y = np.linspace(0.5*Delta, Lt-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

# Rossby number (function of latitude)
Ro = Rom/(1 + Rom*beta*(yc-0.5*Lt))

# choose iteration (-1: last) and layer (0: upper layer)
nt = -1
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
