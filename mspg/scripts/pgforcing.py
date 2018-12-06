#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import def_radius

plt.ion()

file_qg = '../pf_qg.bas0064'
file_dh_qg = '../../msqg/dh.bin'

## qg config

uref = 0.1  # m/s
lref = 50e3 # m

dh_qg = np.fromfile(file_dh_qg,'f4')
h_qg = -np.cumsum(dh_qg)

dh2_qg = 0.5*(dh_qg[1:] + dh_qg[:-1])
h2_qg = -np.cumsum(dh2_qg)
nl_qg = len(dh_qg)



## pg config
nl_pg = 30
N = 64
N1 = N + 1

# PG scales
L = 5000e3  # m
H = 5000    # m
beta = 2.0e-11 # 1/m/s
N2 = 1e-6  #  (1/s**2)


# PG scales
L = 5000e3  # m
H = 5000    # m
beta = 2.0e-11 # 1/m/s
N2 = 1e-6  #  (1/s**2)
Bs = N2*H

Ts = beta*L**3/(N2*H**2)
Us = N2*H**2/(beta*L**2)
Ws = N2*H**3/(beta*L**3)
Thetas = Bs/10/2e-4 # 1/g alpha
Kv = N2*H**4/(beta*L**3)
Kh = N2*H**2/(beta*L)
ys = 0.3


dz_pg = H/nl_pg
dzc = H/nl_pg*np.ones(nl_pg)
zc = 0.5*dzc - np.cumsum(dzc)



# load modes
l2m = np.fromfile('l2m.bin','f4').reshape(nl_pg,nl_pg,N,N)
m2l = np.fromfile('m2l.bin','f4').reshape(nl_pg,nl_pg,N,N)

l2mt = np.fromfile('l2mt.bin','f4').reshape(nl_qg,nl_qg,N,N)
m2lt = np.fromfile('m2lt.bin','f4').reshape(nl_qg,nl_qg,N,N)

# load filter
pf_qg = np.fromfile(file_qg,'f4').reshape(nl_qg,N1,N1).transpose(0,2,1)
pf_qg = pf_qg[:,1:,1:]
# dimesions m^2/s^2/s
pf_qg = pf_qg*uref**2*uref/lref

bf_qg = np.diff(pf_qg,1,0)/dh2_qg.reshape(nl_qg-1,1,1)
bf_qg2 = 0*pf_qg
bf_qg2[0,:,:] = bf_qg[0,:,:]
bf_qg2[-1,:,:] = bf_qg[-1,:,:]
bf_qg2[1:-1,:,:] = 0.5*(bf_qg[:-1,:,:] + bf_qg[1:,:,:])


# p_mod = np.zeros((nl_pg,N,N))
# for ny in range(0,N):
#   for nx in range(0,N):
#     p_mod[:nl_qg,ny,nx] = np.dot(l2mt[:,:,ny,nx],pf_qg[:,ny,nx])

# # remove barotropic pressure
# #p_mod[0,:,:] = 0.

# pf_pg = np.zeros((nl_pg,N,N))
# for ny in range(0,N):
#   for nx in range(0,N):
#     pf_pg[:,ny,nx] = np.dot(m2l[:,:,ny,nx],p_mod[:,ny,nx])


b_mod = np.zeros((nl_pg,N,N))
for ny in range(0,N):
  for nx in range(0,N):
    b_mod[:nl_qg,ny,nx] = np.dot(l2mt[:,:,ny,nx],bf_qg2[:,ny,nx])


bf_pg2 = np.zeros((nl_pg,N,N))
for ny in range(0,N):
  for nx in range(0,N):
    bf_pg2[:,ny,nx] = np.dot(m2l[:,:,ny,nx],b_mod[:,ny,nx])



# bf_pg = np.zeros((nl_pg+1,N,N))
# # buoyancy at cell face
# bf_pg[1:-1,:,:] = np.diff(pf_pg,1,0)/dz_pg
# bf_pg[0,:,:] = bf_pg[1,:,:]
# bf_pg[-1,:,:] = bf_pg[-2,:,:]

# #buoyancy at cell center
# bf_pg = 0.5*(bf_pg[1:,:,:] + bf_pg[:-1,:,:])

# adim
#bf_pg_a = bf_pg*Ts/Bs
bf_pg_a = bf_pg2*Ts/Bs

# Export to basilisk format
bf_pg_o = np.zeros((nl_pg+2,N1,N1))
bf_pg_o[1:-1,1:,1:] = bf_pg_a
bf_pg_o[:,0,0] = N
bf_pg_o = np.transpose(bf_pg_o,(0,2,1))
bf_pg_o.astype('f4').tofile('bf_pg.bas')


psi = bf_qg[5,:,:]*Ts/Bs
vmax = 0.5*np.max(np.abs(psi))

plt.figure()
plt.imshow(psi[::-1,:],vmin=-vmax, vmax=vmax,cmap=plt.cm.bwr)
plt.colorbar()  


nx = 50
ny = 50
plt.figure()
plt.plot(bf_qg[:,ny,nx],h2_qg)
# plt.plot(bf_pg[:,ny,nx],zc)
plt.plot(bf_pg2[:,ny,nx],zc)

