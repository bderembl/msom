#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys,glob,os,re
import scipy.io.netcdf as netcdf


plt.ion()

dir0 = "../outdir_"

if len(sys.argv) > 1:
  dir0 = dir0 + str(format(sys.argv[1])).zfill(4) + '/'

exec(open(dir0 + "params.in").read())

flag_tmp = 0

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

allfilesbf = sorted(glob.glob(dir0 + filebf));
allfilesvd = sorted(glob.glob(dir0 + filevd));
allfilesj1 = sorted(glob.glob(dir0 + filej1));
allfilesj2 = sorted(glob.glob(dir0 + filej2));
allfilesj3 = sorted(glob.glob(dir0 + filej3));
allfilesft = sorted(glob.glob(dir0 + fileft));

nb_files  = len(allfilesp);

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

fileppg = 'psipg_' + str(nl) +'l_N' + str(N) + '.bas'
ppg = np.fromfile(dir0 + fileppg,'f4').reshape(nl,N1,N1).transpose(0,2,1)
ppg = ppg[:,1:,1:]

# create grid
Delta = L0/N

x = np.linspace(0.5*Delta, L0 - 0.5*Delta,N)
xc,yc = np.meshgrid(x,x)

# physical parameters
if (Rom < 0) :
  Ro = 0*yc - Rom
else:
  Ro = Rom/(1 + Rom*beta*(yc-0.5*L0))

# new eq
#ppg = ppg*Ro.reshape(1,N,N)

# backgroud state
#upg = -np.gradient(ppg, axis = 1)/Delta*Ro.reshape(1,N,N)
#vpg = np.gradient(ppg, axis = 2)/Delta*Ro.reshape(1,N,N)
upg = -np.gradient(ppg, axis = 1)/Delta
vpg = np.gradient(ppg, axis = 2)/Delta

kepg = 0.5*(upg**2 + vpg**2)

bpg = np.diff(ppg,1,0)/dhi.reshape(nl-1,1,1)/Ro.reshape(1,N,N)
#pepg = 0.5*bpg**2*Fr[:-1,:,:]**2
pepg = 0.5*bpg**2*Fr[:-1,:,:]**2

kepg_s = (kepg*dh.reshape(nl,1,1)).sum()*Delta*Delta
pepg_s = (pepg*dhi.reshape(nl-1,1,1)).sum()*Delta*Delta
  
i0 = 0
imax = nb_files
#imax = 150

ke_all = np.zeros((imax-i0))
pe_all = np.zeros((imax-i0))
kp_all = np.zeros((imax-i0))

ebf_all = np.zeros((imax-i0))
evd_all = np.zeros((imax-i0))
ej1_all = np.zeros((imax-i0))
ej2_all = np.zeros((imax-i0))
ej3_all = np.zeros((imax-i0))
eft_all = np.zeros((imax-i0))

ke_me = np.zeros((nl,N,N))

nco = 0
for ifi in range(i0,imax):

  p  = np.fromfile(allfilesp[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  q  = np.fromfile(allfilesq[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)

  p  = p [:,1:,1:]
  q  = q [:,1:,1:]

  # new eq
  #p = p*Ro.reshape(1,N,N)


  if len(allfilesbf) > 0:
    ebf = np.fromfile(allfilesbf[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
    evd = np.fromfile(allfilesvd[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
    ej1 = np.fromfile(allfilesj1[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
    ej2 = np.fromfile(allfilesj2[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
    ej3 = np.fromfile(allfilesj3[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
    eft = np.fromfile(allfilesft[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)

    ebf = ebf[:,1:,1:]
    evd = evd[:,1:,1:]
    ej1 = ej1[:,1:,1:]
    ej2 = ej2[:,1:,1:]
    ej3 = ej3[:,1:,1:]
    eft = eft[:,1:,1:]

    ebf_all[nco] = (ebf*dh.reshape(nl,1,1)).sum()*Delta*Delta#/dtout
    evd_all[nco] = (evd*dh.reshape(nl,1,1)).sum()*Delta*Delta#/dtout
    ej1_all[nco] = (ej1*dh.reshape(nl,1,1)).sum()*Delta*Delta#/dtout
    ej2_all[nco] = (ej2*dh.reshape(nl,1,1)).sum()*Delta*Delta#/dtout
    ej3_all[nco] = (ej3*dh.reshape(nl,1,1)).sum()*Delta*Delta#/dtout
    eft_all[nco] = (eft*dh.reshape(nl,1,1)).sum()*Delta*Delta#/dtout

  if ediag == 1:
    pt = p + ppg
  else:
    pt = p

  # not exact near boundaries but very convenient
  #u = -np.gradient(pt, axis = 1)/Delta*Ro.reshape(1,N,N)
  #v = np.gradient(pt, axis = 2)/Delta*Ro.reshape(1,N,N)
  u = -np.gradient(pt, axis = 1)/Delta
  v = np.gradient(pt, axis = 2)/Delta
    
  ke = 0.5*(u**2 + v**2)
  ke_me += ke

#  pt = p + ppg


  b = np.diff(pt,1,0)/dhi.reshape(nl-1,1,1)/Ro.reshape(1,N,N)
  #pe = 0.5*b**2*Fr[:-1,:,:]**2
  pe = 0.5*b**2*Fr[:-1,:,:]**2
  
  ke_all[nco] = (ke*dh.reshape(nl,1,1)).sum()*Delta*Delta
  pe_all[nco] = (pe*dhi.reshape(nl-1,1,1)).sum()*Delta*Delta

  #kp_all[nco] = 0.5*(-p*Ro.reshape(1,N,N)*q*dh.reshape(nl,1,1)).sum()*Delta*Delta
  kp_all[nco] = 0.5*(-p*q*dh.reshape(nl,1,1)).sum()*Delta*Delta
  nco += 1
  
tt = np.linspace(i0,imax-1,imax-i0)*dtout
ke_me = ke_me/(imax-i0)

#pe_all = pe_all - pe_all[0]
#ke_all = ke_all - ke_all[0]
#kp_all = kp_all - kp_all[0]

dedt = np.diff(pe_all + ke_all )/dtout

dedt_sum = ebf_all + evd_all + ej1_all + ej2_all + ej3_all + eft_all

plt.figure()
plt.plot(tt,ke_all, label= 'KE')
plt.plot(tt,pe_all, label= 'PE')
plt.plot(tt,ke_all + pe_all, label= 'PE+Ke')
plt.plot(tt,kp_all, label= 'qp')

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.grid()
plt.xlabel("time")
plt.ylabel("E")
plt.legend()

plt.figure()
plt.plot(tt, ebf_all, label= 'bf')
plt.plot(tt, evd_all, label= 'vd')
plt.plot(tt, ej1_all, label= 'j1')
plt.plot(tt, ej2_all, label= 'j2')
plt.plot(tt, ej3_all, label= 'j3')
plt.plot(tt, eft_all, label= 'ft')

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.grid()
plt.xlabel("time")
plt.ylabel("E")
plt.legend()

plt.figure()
plt.plot(tt[1:], dedt_sum[1:], label= 'sum')
plt.plot(tt[1:],dedt,label='dE/dt')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.grid()
plt.xlabel("time")
plt.ylabel("dE/dt")
plt.legend()

n2 = len(ebf_all)//2
print("Background PE = {0:.1e}".format(pepg_s))
print("Background KE = {0:.1e}".format(kepg_s))
print("ave bf = {0:.1e}, vd ={1:.1e}, j1 = {2:.1e}, j2 = {3:.1e}, j3 = {4:.1e}, ft = {5:.1e}".format(ebf_all[n2:].mean(), evd_all[n2:].mean(), ej1_all[n2:].mean(), ej2_all[n2:].mean(), ej3_all[n2:].mean(), eft_all[n2:].mean()))
print("std bf = {0:.1e}, vd ={1:.1e}, j1 = {2:.1e}, j2 = {3:.1e}, j3 = {4:.1e}, ft = {5:.1e}".format(ebf_all[n2:].std(), evd_all[n2:].std(), ej1_all[n2:].std(), ej2_all[n2:].std(), ej3_all[n2:].std(), eft_all[n2:].std()))

#ke_me = np.where(ke_me>40,40,ke_me)

plt.figure()
plt.contourf(ke_me[0,:,:],50,cmap=plt.cm.hot,extend='both')
plt.colorbar()
