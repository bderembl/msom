#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re
import PowerSpec as ps

plt.ion()

#dir0 = "/home/bderembl/work/basilisk/myrun/msom_run/msqg/run04/outdir/"
dir0 = "../outdir_0017/"

filep = 'po*'
fileq = 'qo*'
filef = 'pf*'

# qg scaling
exec(open(dir0 + "params.in").read())
#Lt = 100
#Rom = 0.025
#beta = 0.5

lref = 1 # km

allfilesp = sorted(glob.glob(dir0 + filep));
allfilesq = sorted(glob.glob(dir0 + fileq));
allfilesf = sorted(glob.glob(dir0 + filef));
nb_files  = len(allfilesp);

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

# grid
Delta = L0/N
x = np.linspace(0.5*Delta, L0-0.5*Delta,N)
y = np.linspace(0.5*Delta, L0-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

Ro = Rom/(1 + Rom*beta*(yc-0.5*L0))

u = np.gradient(p, axis = 1)/Delta
v = np.gradient(p, axis = 2)/Delta
# u = np.gradient(p, axis = 1)/Delta*Ro.reshape(1,N,N)
# v = np.gradient(p, axis = 2)/Delta*Ro.reshape(1,N,N)


ke = 0.5*(u**2 + v**2)

b = np.diff(p,1,0)/dhi.reshape(nl-1,1,1)
pe = 0.5*b**2*Fr[:-1,:,:]**2/(Ro.reshape(1,N,N))**2


plt.figure()
for keplot in range(0,2):
  for il in range(0,nl):

    plt.clf()

    if keplot == 1:
      psi = ke[il,:,:]
    else:
      if il < nl-1:
        psi = pe[il,:,:]
      else:
        psi = 0

    #psi = np.sin(2*np.pi*0.1*xc)*np.sin(2*np.pi*0.1*yc)
    if keplot == 0 and il == nl-1:
      psi = 0
      # don't plot
    else:
      kspec,pspec = ps.get_spectrum(psi,xc,yc,window=None,detrend=None)
      lam = 2*np.pi/kspec*lref
      
      
      ki = 1.
      k0 = np.argmin(np.abs(kspec-ki))
      ps0 = pspec[k0]
      
      k3 = np.array([3e-1,3e0])
      s3 = (k3/ki)**-(3)*ps0 # k-5  QG
      lam3 = 2*np.pi/k3*lref
      
      k5 = np.array([3e-1,3e0])
      s5 = (k5/ki)**-(5/3)*ps0 # k-5  QG
      lam5 = 2*np.pi/k5*lref
      
      #plt.loglog(lam,pspec,'k-',label=r'PE')
      #plt.loglog(lam3,s3,'k--',label=r'$k^{-3}$')
      #plt.loglog(lam5,s5,'k-.',label=r'$k^{-5/3}$')
      if keplot == 1:
        tag = 'KE'+str(il+1)
      else:
        tag = 'PE'+str(il+1)
      plt.loglog(kspec,pspec,'-',label=tag)
#      plt.loglog(kspec,pspec,'r-',label=tag)
      
      plt.loglog(k3,s3,'k--',label=r'$k^{-3}$', linewidth=1)
      plt.loglog(k5,s5,'k-.',label=r'$k^{-5/3}$', linewidth=1)
      plt.xlabel('k')
      plt.legend()
      plt.grid()
      print(tag)

      plt.savefig('figures/' + tag + '.png', bbox_inches='tight')
  
#plt.gca().invert_xaxis()
#plt.savefig('figures/' + 'PE_all' +'.pdf', bbox_inches='tight')
