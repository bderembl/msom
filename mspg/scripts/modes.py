#!/usr/bin/env python

import glob,os,re
import numpy as np
import matplotlib.pyplot as plt
import def_radius
import spoisson

plt.ion()

dir0 = "../outdir/"

file_b = 'b*'
file_u = 'u*'

flag_uniform_strat = 0 # 0: average froude number, 1: variable
adjust_psi_coef = 0.1  # multiply psi by coef to study weaker flow

allfilesu = sorted(glob.glob(dir0 + file_u));
allfilesb = sorted(glob.glob(dir0 + file_b));
nb_files  = len(allfilesu);

# dimensions
b = np.fromfile(allfilesb[0],'f4')
N = int(b[0])
N1 = N + 1
nl2 = int(len(b)/N1**2)
nl = nl2 - 2


# PG scales
L = 5000e3  # m
H = 5000    # m
beta = 2.0e-11 # 1/m/s
N2 = 1e-6  #  (1/s**2)

Ts = beta*L**3/(N2*H**2)
Us = N2*H**2/(beta*L**2)
Ws = N2*H**3/(beta*L**3)
Bs = N2*H
Thetas = Bs/10/2e-4 # 1/g alpha
Kv = N2*H**4/(beta*L**3)
Kh = N2*H**2/(beta*L)
ys = 0.3

# qg scales
u_qg = 0.1  # m/s
l_qg = 50e3 # m

# grid
Delta = 1/N
Deltad = L*Delta
x = np.linspace(0.5*Delta, 1-0.5*Delta,N)
y = ys + np.linspace(0.5*Delta, 1-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

#temporary
dz = H/nl*np.ones(nl)
z = 0.5*dz - np.cumsum(dz)

# compute mean
b = np.zeros((nl2,N1,N1))
uv = np.zeros((2*nl2,N1,N1))

nme = 10
for ifi in range(nb_files-nme,nb_files): 
  b += np.fromfile(allfilesb[-1],'f4').reshape(nl2,N1,N1).transpose(0,2,1)
  uv += np.fromfile(allfilesu[-1],'f4').reshape(2*nl2,N1,N1).transpose(0,2,1)

b = b[1:-1,1:,1:]/nme
u = uv[2:-2:2,1:,1:]/nme
v = uv[3:-2:2,1:,1:]/nme


# adjust stratification above threshold
N2_min = 2e-7
N2_l = Bs*np.diff(b,1,0)/np.diff(z).reshape(nl-1,1,1)
for nz in range(0,nl-1):
  for nx in range(0,N):
    for ny in range(0,N):
      if (N2_l[nz,ny,nx] <  N2_min):
        dN = N2_min - N2_l[nz,ny,nx]
        b[nz+1:,ny,nx] -= (dN*(z[nz] - z[nz+1]))/Bs

#gpmin = 3e-4
gp = -Bs*np.diff(b,1,0)
#gp = np.where(gp<gpmin,gpmin,gp)
f0 = yc*L*beta

# il: layer interface
il = [0,3,8,16,30]        # 4 layers
#il = [0,3,6,9,12,15,18,21,24,27,30] # 10 layers
#il = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]  # 15 layers
#il = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,28,30] # 23 layers
#il = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30] # 30 layers
nlt = len(il)-1
bt = np.zeros((nlt,N,N))
ut = np.zeros((nlt,N,N))
vt = np.zeros((nlt,N,N))
gpt = np.zeros((nlt-1,N,N))
dzt = np.zeros((nlt))

for n2 in range(0,nlt):
  bt[n2,:,:] = np.mean(b[il[n2]:il[n2+1],:,:],0)
  ut[n2,:,:] = np.mean(u[il[n2]:il[n2+1],:,:],0)
  vt[n2,:,:] = np.mean(v[il[n2]:il[n2+1],:,:],0)
gpt[:,:,:] = -Bs*np.diff(bt,1,0)
#gpt = np.where(gpt<gpmin,gpmin,gpt)

if flag_uniform_strat:
  gpt = 0*gpt + gpt.mean(1).mean(1).reshape((nlt-1,1,1))

for n2 in range(0,nlt):
  dzt[n2] = np.sum(dz[il[n2]:il[n2+1]])

dzi = 0.5*(dzt[:-1] + dzt[1:])
zt = 0.5*dzt - np.cumsum(dzt)

N2lt = gpt/dzi.reshape((nlt-1,1,1))

rd = np.zeros((nl,N,N))
l2m = np.zeros((nl,nl,N,N))
m2l = np.zeros((nl,nl,N,N))

for nx in range(0,N):
  for ny in range(0,N):
    rd[:,ny,nx] = def_radius.cal_rad(dz,gp[:,ny,nx],f0[ny,nx])
    l2,m2 = def_radius.cal_transfo(dz,gp[:,ny,nx],f0[ny,nx])
    l2m[:,:,ny,nx] = l2
    m2l[:,:,ny,nx] = m2

l2m.astype('f4').tofile('l2m.bin')
m2l.astype('f4').tofile('m2l.bin')

rdt = np.zeros((nlt,N,N))
l2mt = np.zeros((nlt,nlt,N,N))
m2lt = np.zeros((nlt,nlt,N,N))

for nx in range(0,N):
  for ny in range(0,N):
    rdt[:,ny,nx] = def_radius.cal_rad(dzt,gpt[:,ny,nx],f0[ny,nx])
    l2,m2 = def_radius.cal_transfo(dzt,gpt[:,ny,nx],f0[ny,nx])
    l2mt[:,:,ny,nx] = l2
    m2lt[:,:,ny,nx] = m2

# # sort
# rdt_s = 0*rdt
# l2mt_s = 0*l2mt
# m2lt_s = 0*m2lt
# for nx in range(0,N):
#   for ny in range(0,N):


l2mt.astype('f4').tofile('l2mt.bin')
m2lt.astype('f4').tofile('m2lt.bin')

psi_ls = np.zeros((nlt,N1,N1))
# compute large scale stream function
for n2 in range(0,nlt):
  #fu = f0*Us*ut[n2,:,:] 
  #fv = f0*Us*vt[n2,:,:] 
  fu = Us*ut[n2,:,:] 
  fv = Us*vt[n2,:,:] 
  zeta = (fv[1:,1:] - fv[1:,:-1] - fu[1:,1:] + fu[:-1,1:])/Deltad
  psi_ls[n2,1:-1,1:-1] = Deltad**2*spoisson.sol(zeta)

ir = 1

ci = [1,2,5,10,20,30,40,50,60,70,80,90]
plt.figure()
plt.contourf(x,y,b[0,:,:])
plt.colorbar()
CS = plt.contour(x,y,rd[ir,:,:]*1e-3, ci, colors='k')
plt.clabel(CS, inline=1, fontsize=10)

plt.figure()
plt.contourf(x,y,bt[0,:,:])
plt.colorbar()
CS = plt.contour(x,y,rdt[ir,:,:]*1e-3, ci, colors='k')
plt.clabel(CS, inline=1, fontsize=10)

nx = 50
ny = 50

plt.figure()
for nm in range(0,nlt):
  plt.plot(m2l[:,nm,ny,nx],z)

plt.figure()
for nm in range(0,nlt):
  plt.plot(m2lt[:,nm,ny,nx],zt)


u_mod = np.zeros(nl)
u_mod[:nlt] = np.dot(l2mt[:,:,ny,nx],ut[:,ny,nx])
u_proj = np.dot(m2l[:,:,ny,nx],u_mod)

u_mod2 = np.dot(l2m[:,:,ny,nx],u[:,ny,nx])
u_mod2[nlt:] = 0.
u_proj2 = np.dot(m2l[:,:,ny,nx],u_mod2)

plt.figure()
plt.plot(u[:,ny,nx],z,'k+-')
plt.plot(ut[:,ny,nx],zt,'k+--')
plt.plot(u_proj,z,'r')
plt.plot(u_proj2,z,'b')

b_mod = np.zeros(nl)
b_mod[:nlt] = np.dot(l2mt[:,:,ny,nx],bt[:,ny,nx])
b_proj = np.dot(m2l[:,:,ny,nx],b_mod)

b_mod2 = np.dot(l2m[:,:,ny,nx],b[:,ny,nx])
b_mod2[nlt:] = 0.
b_proj2 = np.dot(m2l[:,:,ny,nx],b_mod2)

plt.plot(b[:,ny,nx],z,'k+-')
plt.plot(bt[:,ny,nx],zt,'k+--')
plt.plot(b_proj,z,'r')
plt.plot(b_proj2,z,'b')

# adjust psi
psi_ls = adjust_psi_coef*psi_ls

# adim
psi_ls_a = psi_ls/(l_qg*u_qg)
gpt_a    = gpt*l_qg/u_qg**2
dzt_a    = dzt/l_qg
Fr = u_qg/(np.sqrt(N2lt)*H)
dht_a    = dzt/H

# Export to basilisk format
psi_ls_o = 0*psi_ls
psi_ls_o[:,1:,1:] = 0.25*(psi_ls_a[:,:-1,:-1] + psi_ls_a[:,1:,1:] + psi_ls_a[:,:-1,1:] + psi_ls_a[:,1:,:-1])
psi_ls_o[:,0,:] = 0
psi_ls_o[:,:,0] = 0
psi_ls_o[:,0,0] = N
psi_ls_o = np.transpose(psi_ls_o,(0,2,1))

fileppg = 'psipg_' + str(nlt) +'l.bas'
psi_ls_o.astype('f4').tofile(fileppg)

Fr_o = np.zeros((nlt,N1,N1))
Fr_o[:-1,1:,1:] = Fr
Fr_o[:,0,:] = 0
Fr_o[:,:,0] = 0
Fr_o[:,0,0] = N
Fr_o = np.transpose(Fr_o,(0,2,1))
fileFr = 'frpg_' + str(nlt) +'l.bas'
Fr_o.astype('f4').tofile(fileFr)

# gp_o = np.zeros((nlt,N1,N1))
# gp_o[:-1,1:,1:] = gpt_a
# gp_o[:,0,:] = 0
# gp_o[:,:,0] = 0
# gp_o[:,0,0] = N
# gp_o = np.transpose(gp_o,(0,2,1))
# gp_o.astype('f4').tofile('gppg.bas')


plt.figure()
plt.contourf(xc,yc,psi_ls_o[0,1:,1:].T,20)
plt.colorbar()

fileh = 'dh_' + str(nlt) +'l.bin'
#dzt_a.astype('f4').tofile('dh_old.bin')
dht_a.astype('f4').tofile(fileh)

plt.figure()
plt.contourf(xc,yc,ut[0,:,:],20)
plt.colorbar()

print('Rom = {0:.2e}'.format(u_qg/f0.mean()/l_qg))
