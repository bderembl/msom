#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pickle

plt.ion()

file0 = 'params.dat'
file1 = 'sol.dat'
file2 = 'cont_diag.dat'

with open(file0,'rb') as fout:
  params = pickle.loads(fout.read())


data = np.loadtxt(file1)
cont_par = np.loadtxt(file2)

# N and nl have to match the dimensions in pg.c
N = params["N"]
nl = params["nl"]
N1 = N + 1

s = np.linspace(0, -1, nl)

#si_t,ne = data.shape


# choose vertical level (0=top, nl-1=bottom)
l = 0

#plot velocity vector every nsk pts
nsk = 5

nt = -1

b = data[nt,:N*N*nl].reshape((N,N,nl)).T
u = data[nt,N*N*nl:N*N*nl+N1*N*nl].reshape((N1,N,nl)).T
v = data[nt,N*N*nl+N1*N*nl:].reshape((N,N1,nl)).T

u = u[:,:,:-1]
v = v[:,:-1,:]

vmin = np.min(b)
vmax = np.max(b)

plt.figure()
plt.contourf(b[0,:,:])
#plt.contourf(b[:,:,15])  
plt.colorbar()


plt.figure()
plt.contourf(b[1,:,:])
#plt.contourf(b[:,:,15])  
plt.colorbar()


plt.figure()
plt.contourf(b[:,:,15])  
plt.colorbar()
