#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re,sys
import fftlib as myfftlib

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
nl = int(len(b)/N1**2)

fileh = 'dh_' + str(nl) +'l.bin'

dh = np.fromfile(dir0 + fileh,'f4')
dhi = 0.5*(dh[:-1] + dh[1:])

fileFr = 'frpg_' + str(nl) +'l_N' + str(N) + '.bas'

nt = -1
il = 0

p  = np.fromfile(allfilesp[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
q  = np.fromfile(allfilesq[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
pf = np.fromfile(allfilesf[nt],'f4').reshape(nl,N1,N1).transpose(0,2,1)
Fr = np.fromfile(dir0 + fileFr,'f4').reshape(nl,N1,N1).transpose(0,2,1)

p  = p[:,1:,1:]
q  = q[:,1:,1:]
pf = pf[:,1:,1:]
Fr = Fr[:,1:,1:]

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

u = np.gradient(p, axis = 1)/Delta
v = np.gradient(p, axis = 2)/Delta
# u = np.gradient(p, axis = 1)/Delta*Ro.reshape(1,N,N)
# v = np.gradient(p, axis = 2)/Delta*Ro.reshape(1,N,N)

omega = np.gradient(np.gradient(p, axis = 1), axis = 1)/Delta**2 \
  +     np.gradient(np.gradient(p, axis = 2), axis = 2)/Delta**2

ke = 0.5*(u**2 + v**2)

b = np.diff(p,1,0)/dhi.reshape(nl-1,1,1)
pe = 0.5*b**2*Fr[:-1,:,:]**2/(Ro.reshape(1,N,N))**2

nmax = nl # 1: just surface level, nl: all levels
flag_plot = [2,3,4,5,6,7] # 0: PE, 1: KE, 2: 
flag_plot = [1] # 0: PE, 1: KE, 2: 

Nkr = myfftlib.get_len_wavenumber(N, Delta)

plt.figure()
for plotvar in flag_plot:
  eflux = np.zeros((nl,Nkr))
  for il in range(0,nmax):

    #plt.clf()

    if plotvar == 0: # PE
      if il < nl-1:
        psi = pe[il,:,:]
      else:
        psi = 0
    elif plotvar == 1: # KE
#      psi = ke[il,:,:]
      psi = p[il,:,:]*omega[il,:,:]
    elif plotvar == 2: # bf
      psi = ebf[il,:,:]
    elif plotvar == 3: # vd
      psi = evd[il,:,:]
    elif plotvar == 4: # j1
      psi = ej1[il,:,:]
    elif plotvar == 5: # j2
      psi = ej2[il,:,:]
    elif plotvar == 6: # j3
      psi = ej3[il,:,:]
    elif plotvar == 7: # ft
      psi = eft[il,:,:]

    #psi = np.sin(2*np.pi*0.1*xc)*np.sin(2*np.pi*0.1*yc)
    if plotvar == 0 and il == nl-1:
      psi = 0
      # don't plot
    else:
      kspec,pspec = myfftlib.get_spec_1D(psi,psi,Delta)
      kspec,flux = myfftlib.get_flux(-p[il,:,:],psi,Delta)
      eflux[il,:] = flux*dh[il]      
      
      ki = 1.
      k0 = np.argmin(np.abs(kspec-ki))
      ps0 = pspec[k0]
      
      k3 = np.array([3e-1,3e0])
      s3 = (k3/ki)**-(3)*ps0 # k-5  QG
      
      k5 = np.array([3e-1,3e0])
      s5 = (k5/ki)**-(5/3)*ps0 # k-5  QG
      
      if plotvar == 0:
        tag = 'PE'+str(il+1)
      elif plotvar == 1:
        tag = 'KE'+str(il+1)
      elif plotvar == 2: # bf
        tag = 'bf'+str(il+1)
      elif plotvar == 3: # vd
        tag = 'vd'+str(il+1)
      elif plotvar == 4: # j1
        tag = 'j1'+str(il+1)
      elif plotvar == 5: # j2
        tag = 'j2'+str(il+1)
      elif plotvar == 6: # j3
        tag = 'j3'+str(il+1)
      elif plotvar == 7: # ft
        tag = 'ft'+str(il+1)

#      plt.semilogx(kspec,flux,'-',label=tag)
      
      plt.loglog(kspec,pspec,'r-',label=tag)      
      plt.loglog(k3,s3,'k--',label=r'$k^{-3}$', linewidth=1)
      plt.loglog(k5,s5,'k-.',label=r'$k^{-5/3}$', linewidth=1)
      plt.xlabel('k (cycle/l)')
      plt.legend()
      plt.grid()
      print(tag)
      plt.savefig('figures/' + tag + '.png', bbox_inches='tight')
  
#plt.semilogx(kspec,np.sum(eflux,axis=0),'-', label='sum');
plt.legend()

#plt.gca().invert_xaxis()
#plt.savefig('figures/' + 'PE_all' +'.pdf', bbox_inches='tight')
