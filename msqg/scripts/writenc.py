#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys,glob,os,re
import scipy.io.netcdf as netcdf

dtflt = -1

#plt.ion()


dir0 = "../outdir_"

if len(sys.argv) > 1:
  dir0 = dir0 + str(format(sys.argv[1])).zfill(4) + '/'

exec(open(dir0 + "params.in").read())

filep = 'po*'
fileq = 'qo*'
filef = 'pf*'
filet = 'ptr*'

filebf = 'de_bf*'
filevd = 'de_vd*'
filej1 = 'de_j1*'
filej2 = 'de_j2*'
filej3 = 'de_j3*'
fileft = 'de_ft*'

allfilesp = sorted(glob.glob(dir0 + filep));
allfilesq = sorted(glob.glob(dir0 + fileq));
allfilesf = sorted(glob.glob(dir0 + filef));

allfilest = sorted(glob.glob(dir0 + filet));

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

# create netcdf file
fileout = dir0 + 'vars.nc'
f = netcdf.netcdf_file(fileout,'w')

f.createDimension('t',None)
f.createDimension('z',nl)
f.createDimension('y',N)
f.createDimension('x',N)

tpo = f.createVariable('t', 'f', ('t',))
zpo = f.createVariable('z', 'f', ('z',))
ypo = f.createVariable('y', 'f', ('y',))
xpo = f.createVariable('x', 'f', ('x',))

po  = f.createVariable('p' , 'f', ('t','z','y','x',))
qo  = f.createVariable('q' , 'f', ('t','z','y','x',))
if len(allfilesf) > 0:
  pof  = f.createVariable('pf' , 'f', ('t','z','y','x',))
if len(allfilesbf) > 0:
  ebfo = f.createVariable('ebf' , 'f', ('t','z','y','x',))
  evdo = f.createVariable('evd' , 'f', ('t','z','y','x',))
  ej1o = f.createVariable('ej1' , 'f', ('t','z','y','x',))
  ej2o = f.createVariable('ej2' , 'f', ('t','z','y','x',))
  ej3o = f.createVariable('ej3' , 'f', ('t','z','y','x',))
  efto = f.createVariable('eft' , 'f', ('t','z','y','x',))
if len(allfilest) > 0:
  ptro = f.createVariable('ptr' , 'f', ('t','z','y','x',))


zpo[:] = np.arange(nl)
ypo[:] = np.arange(N)
xpo[:] = np.arange(N)

# tpo[0] = 1.0
# bo[0,:,:,:] = 1.0

ifi = nb_files - 1

for ifi in range(0,nb_files):
#for ifi in range(0,100):

  p  = np.fromfile(allfilesp[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  q  = np.fromfile(allfilesq[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)

  p  = p [:,1:,1:]
  q  = q [:,1:,1:]

  po [ifi,:,:,:] = p [:,:,:]
  qo [ifi,:,:,:] = q [:,:,:]

  if len(allfilesf) > 0:
    pf = np.fromfile(allfilesf[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
    pf = pf[:,1:,1:]
    pof[ifi,:,:,:] = pf[:,:,:]

  if len(allfilest) > 0:
    ptr = np.fromfile(allfilest[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
    ptr = ptr[:,1:,1:]
    ptro[ifi,:,:,:] = ptr[:,:,:]


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

    ebfo[ifi,:,:,:] = ebf[:,:,:]
    evdo[ifi,:,:,:] = evd[:,:,:]
    ej1o[ifi,:,:,:] = ej1[:,:,:]
    ej2o[ifi,:,:,:] = ej2[:,:,:]
    ej3o[ifi,:,:,:] = ej3[:,:,:]
    efto[ifi,:,:,:] = eft[:,:,:]

  tpo[ifi] = dtout*ifi

f.close()
print ("nb points in file: {0}".format(nb_files))

