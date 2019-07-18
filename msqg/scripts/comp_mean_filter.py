#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re
import scipy.io.netcdf as netcdf


plt.ion()

dir0 = "../outdir_0014/"
#dir0 = "outdir_filt15_newpg/"
exec(open(dir0 + "params.in").read())

filep = 'po*'
fileq = 'qo*'
filef = 'pf*'

allfilesp = sorted(glob.glob(dir0 + filep));
allfilesq = sorted(glob.glob(dir0 + fileq));
allfilesf = sorted(glob.glob(dir0 + filef));
nb_files  = len(allfilesp);

b = np.fromfile(allfilesp[0],'f4')
N = int(b[0])
N1 = N + 1
nl = int(len(b)/N1**2)

x = L0*np.arange(N)/N


fileppg = 'psipg_' + str(nl) +'l_N' + str(N) + '.bas'
ppg = np.fromfile(dir0 + fileppg,'f4').reshape(nl,N1,N1).transpose(0,2,1)
ppg = ppg[:,1:,1:]


p_me  = np.zeros((nl,N,N))
pf_me = np.zeros((nl,N,N))
nme = 0
ifi = nb_files - 1
for ifi in range(0,nb_files):
#for ifi in range(0,100):

  p  = np.fromfile(allfilesp[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  pf = np.fromfile(allfilesf[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)

  p  = p[:,1:,1:]
  pf = pf[:,1:,1:]

  p_me  += p
  pf_me += pf
  nme += 1


p_me  /= nme
pf_me /= nme

pf_o = np.transpose(pf_me,(0,2,1))
pf_o.astype('f4').tofile('pf_qg.bas')

psi = ppg[0,1:,1:]
#psi = ppg[0,1:,1:] + pf_me[0,1:,1:]
vmax = 0.5*np.max(np.abs(psi))

plt.figure()
plt.contour(psi,15, colors='k', linewidths=1)
