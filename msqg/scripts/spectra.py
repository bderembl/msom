#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re,sys
import fftlib as myfftlib
from matplotlib import colors

plt.ion()

#dir0 = "/home/bderembl/work/basilisk/myrun/msom_run/msqg/run04/outdir/"
dir0 = "../outdir_"

if len(sys.argv) > 1:
  dir0 = dir0 + str(format(sys.argv[1])).zfill(4) + '/'


filep = 'po*'
fileq = 'qo*'
filef = 'pf*'

filebf = 'de_bf*'
filevd = 'de_vd*'
filej1 = 'de_j1*'
filej2 = 'de_j2*'
filej3 = 'de_j3*'
fileft = 'de_ft*'

# qg scaling
exec(open(dir0 + "params.in").read())

lref = 1 # km

allfilesp = sorted(glob.glob(dir0 + filep));
allfilesq = sorted(glob.glob(dir0 + fileq));
allfilesf = sorted(glob.glob(dir0 + filef));
nb_files  = len(allfilesp);

allfilesbf = sorted(glob.glob(dir0 + filebf));
allfilesvd = sorted(glob.glob(dir0 + filevd));
allfilesj1 = sorted(glob.glob(dir0 + filej1));
allfilesj2 = sorted(glob.glob(dir0 + filej2));
allfilesj3 = sorted(glob.glob(dir0 + filej3));
allfilesft = sorted(glob.glob(dir0 + fileft));

b = np.fromfile(allfilesp[0],'f4')
N = int(b[0])
N1 = N + 1
N2 = N + 2
nl = int(len(b)/N1**2)

fileh = 'dh_' + str(nl) +'l.bin'

dh = np.fromfile(dir0 + fileh,'f4')
dhi = 0.5*(dh[:-1] + dh[1:])

fileFr = 'frpg_' + str(nl) +'l_N' + str(N) + '.bas'

nt = -10
il = 0

p  = np.fromfile(allfilesp[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
q  = np.fromfile(allfilesq[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
Fr = np.fromfile(dir0 + fileFr,'f4').reshape(nl,N1,N1).transpose(0,2,1)

p  = p[:,1:,1:]
q  = q[:,1:,1:]
Fr = Fr[:,1:,1:]

# padded p
po2 = np.zeros((nl,N2,N2))
po2[:,1:-1,1:-1] = p

po2[:,0,:]  = -po2[:,1,:]
po2[:,-1,:] = -po2[:,-2,:]
po2[:,:,0]  = -po2[:,:,1]
po2[:,:,-1] = -po2[:,:,-2]

# corners
po2[:,0,0]   = -po2[:,0,1]   - po2[:,1,0]   - po2[:,1,1]
po2[:,-1,0]  = -po2[:,-1,1]  - po2[:,-2,0]  - po2[:,-2,1]
po2[:,0,-1]  = -po2[:,1,-1]  - po2[:,0,-2]  - po2[:,1,-2]
po2[:,-1,-1] = -po2[:,-1,-2] - po2[:,-2,-2] - po2[:,-2,-1]



if len(allfilesbf) > 0:
  ebf = np.fromfile(allfilesbf[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  evd = np.fromfile(allfilesvd[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  ej1 = np.fromfile(allfilesj1[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  ej2 = np.fromfile(allfilesj2[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  ej3 = np.fromfile(allfilesj3[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  eft = np.fromfile(allfilesft[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  
  ebf = ebf[:,1:,1:]
  evd = evd[:,1:,1:]
  ej1 = ej1[:,1:,1:]
  ej2 = ej2[:,1:,1:]
  ej3 = ej3[:,1:,1:]
  eft = eft[:,1:,1:]  

# grid
Delta = L0/N
x = np.linspace(0.5*Delta, L0-0.5*Delta,N)
y = np.linspace(0.5*Delta, L0-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

if (Rom < 0) :
  Ro = 0*yc - Rom
else:
  Ro = Rom/(1 + Rom*beta*(yc-0.5*L0))

u = -(po2[:,2:,1:-1] - po2[:,:-2,1:-1])/(2*Delta)
v =  (po2[:,1:-1,2:] - po2[:,1:-1,:-2])/(2*Delta)

omega = (po2[:,2:,1:-1] + po2[:,:-2,1:-1] + po2[:,1:-1,2:] + po2[:,1:-1,:-2] - 4*p)/Delta**2

ke = 0.5*(u**2 + v**2)

b = np.diff(p,1,0)/dhi.reshape(nl-1,1,1)/Ro.reshape(1,N,N)
pe = 0.5*b**2*Fr[:-1,:,:]**2

Nkr = myfftlib.get_len_wavenumber(N, Delta)
k,l,K,kr = myfftlib.get_wavenumber(N,Delta)

kespec = np.zeros((nl,Nkr))
pespec = np.zeros((nl,Nkr))

for il in range(0,nl):

#  kepspec = myfftlib.get_spec_1D(-p[il,:,:],omega[il,:,:],Delta)
#  kespec[il,:] = 0.5*kepspec*dh[il]      
  upspec = myfftlib.get_spec_1D(u[il,:,:],u[il,:,:],Delta)
  vpspec = myfftlib.get_spec_1D(v[il,:,:],v[il,:,:],Delta)
  kespec[il,:] = 0.5*(upspec + vpspec)*dh[il]      

  if il < nl-1:
    pepspec = myfftlib.get_spec_1D(b[il,:,:]*Fr[il,:,:],b[il,:,:]*Fr[il,:,:],Delta)
    pespec[il,:] = 0.5*pepspec*dhi[il]
      
ke_all = np.sum(kespec,axis=0)
pe_all = np.sum(pespec,axis=0)

ki = 1.
k0 = np.argmin(np.abs(kr-ki))
ps0 = ke_all[k0]

k5 = np.array([3e-1,3e0])
s5 = (k5/ki)**-(5/3)*ps0 # k-5/3  SQG

ki = 0.2
k0 = np.argmin(np.abs(kr-ki))
ps0 = ke_all[k0]

k3 = np.array([1e-1,1.])
s3 = (k3/ki)**-(3)*ps0 # k-3  QG


plt.figure()
plt.loglog(kr,ke_all,'k-',label="KE")
#plt.loglog(kr,kespec[0,:],'k-',label="KE")
plt.loglog(kr,pe_all,'r-',label="PE")
plt.loglog(k3,s3,'k--',label=r'$k^{-3}$', linewidth=1)
plt.loglog(k5,s5,'k-.',label=r'$k^{-5/3}$', linewidth=1)
plt.xlabel('|K| (cycle/l)')
plt.ylabel('E density (l^3/t^2)')
plt.legend()
plt.grid()

#k = np.delete(k,int(N/2),axis=0)
#k = np.delete(k,int(N/2),axis=1)
#l = np.delete(l,int(N/2),axis=0)
#l = np.delete(l,int(N/2),axis=1)
#
#klog = np.sign(k)*np.log10(np.abs(k))
#llog = np.sign(l)*np.log10(np.abs(l))
#
#il = 0
#spec2d_u = myfftlib.get_spec_2D(u[il,:,:],u[il,:,:],Delta)
#spec2d_u = np.delete(spec2d_u,int(N/2),axis=1)
#spec2d_u = np.delete(spec2d_u,int(N/2),axis=0)
#spec2d_u = np.where(spec2d_u<1,1,spec2d_u)
#
#spec2d_v = myfftlib.get_spec_2D(v[il,:,:],v[il,:,:],Delta)
#spec2d_v = np.delete(spec2d_v,int(N/2),axis=1)
#spec2d_v = np.delete(spec2d_v,int(N/2),axis=0)
#spec2d_v = np.where(spec2d_v<1,1,spec2d_v)
#
#spec2d_ke = 0.5*(spec2d_u + spec2d_v)*dhi[il]

#plt.figure()
#plt.contourf(k,l,spec2d_u + spec2d_v,norm = colors.LogNorm())
#plt.xscale('symlog')
#plt.yscale('symlog')
#plt.colorbar()

#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
