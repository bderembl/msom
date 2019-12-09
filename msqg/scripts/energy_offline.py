#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re,sys
sys.path.append('../')
import qg as bas
import fftlib as myfftlib


dir0 = "./"
dir0 = "../outdir_"

if len(sys.argv) > 1:
  dir0 = dir0 + str(format(sys.argv[1])).zfill(4) + '/'

# read parameters
exec(open(dir0 + "params.in").read())

# needed for basilisk
dirtmp = os.getcwd()
os.chdir(dir0)

bas.read_params()
bas.init_grid(N)
bas.set_vars()
bas.set_vars_energy()
bas.set_const()

os.chdir(dirtmp)

p  = np.zeros((nl,N,N))
q  = np.zeros((nl,N,N))
bf = np.zeros((nl,N,N))
vd = np.zeros((nl,N,N))
j1 = np.zeros((nl,N,N))
j2 = np.zeros((nl,N,N))
j3 = np.zeros((nl,N,N))
ft = np.zeros((nl,N,N))

filep = 'po*'
fileq = 'qo*'
filef = 'pf*'

filebf = 'de_bf*'
filevd = 'de_vd*'
filej1 = 'de_j1*'
filej2 = 'de_j2*'
filej3 = 'de_j3*'
fileft = 'de_ft*'

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
nl = int(len(b)/N1**2)

fileh = 'dh_' + str(nl) +'l.bin'

dh = np.fromfile(dir0 + fileh,'f4')
dhi = 0.5*(dh[:-1] + dh[1:])

fileFr = 'frpg_' + str(nl) +'l_N' + str(N) + '.bas'
Fr = np.fromfile(dir0 + fileFr,'f4').reshape(nl,N1,N1).transpose(0,2,1)
Fr = Fr[:,1:,1:]

# grid
Delta = L0/N
x = np.linspace(0.5*Delta, L0-0.5*Delta,N)
y = np.linspace(0.5*Delta, L0-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

if (Rom < 0) :
  Ro = 0*yc - Rom
else:
  Ro = Rom/(1 + Rom*beta*(yc-0.5*L0))


Nkr = myfftlib.get_len_wavenumber(N,Delta)

eflux_bf = np.zeros((nl,Nkr))
eflux_vd = np.zeros((nl,Nkr))
eflux_j1 = np.zeros((nl,Nkr))
eflux_j2 = np.zeros((nl,Nkr))
eflux_j3 = np.zeros((nl,Nkr))
eflux_ft = np.zeros((nl,Nkr))

nme = 0

nt0 = -2
nt1 = -1
for nt in range (nt0,nt1):

  p  = np.fromfile(allfilesp[-1],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  q  = np.fromfile(allfilesq[-1],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  p  = p[:,1:,1:]
  q  = q[:,1:,1:]

  bas.pystep(p,bf,vd,j1,j2,j3,ft,1)
  print("end step")
  nme += 1

  for il in range(0,nl):
    print("layer {0}".format(il))
    kspec,flux_bf  = myfftlib.get_flux(-p[il,:,:],bf[il,:,:],Delta)
    kspec,flux_vd  = myfftlib.get_flux(-p[il,:,:],vd[il,:,:],Delta)
    kspec,flux_j1  = myfftlib.get_flux(-p[il,:,:],j1[il,:,:],Delta)
    kspec,flux_j2  = myfftlib.get_flux(-p[il,:,:],j2[il,:,:],Delta)
    kspec,flux_j3  = myfftlib.get_flux(-p[il,:,:],j3[il,:,:],Delta)
    kspec,flux_ft  = myfftlib.get_flux(-p[il,:,:],ft[il,:,:],Delta)
    
    eflux_bf[il,:] +=  flux_bf*dh[il]
    eflux_vd[il,:] +=  flux_vd*dh[il]
    eflux_j1[il,:] +=  flux_j1*dh[il]
    eflux_j2[il,:] +=  flux_j2*dh[il]
    eflux_j3[il,:] +=  flux_j3*dh[il]
    eflux_ft[il,:] +=  flux_ft*dh[il]

eflux_bf /= nme
eflux_vd /= nme
eflux_j1 /= nme
eflux_j2 /= nme
eflux_j3 /= nme
eflux_ft /= nme

plt.figure()
plt.plot(np.log10(kspec), np.sum(eflux_bf,0), label= 'bf')
plt.plot(np.log10(kspec), np.sum(eflux_vd,0), label= 'vd')
plt.plot(np.log10(kspec), np.sum(eflux_j1,0), label= 'j1')
plt.plot(np.log10(kspec), np.sum(eflux_j2,0), label= 'j2')
plt.plot(np.log10(kspec), np.sum(eflux_j3,0), label= 'j3')
plt.plot(np.log10(kspec), np.sum(eflux_ft,0), label= 'ft')

plt.legend()
plt.show()

bas.trash_vars()
bas.trash_vars_energy()
