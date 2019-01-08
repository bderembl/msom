#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob,os,re
import scipy.io.netcdf as netcdf


#plt.ion()

dir0 = "../outdir/"

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
pof  = f.createVariable('pf' , 'f', ('t','z','y','x',))

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
  pf = np.fromfile(allfilesf[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)

  p  = p [:,1:,1:]
  q  = q [:,1:,1:]
  pf = pf[:,1:,1:]
  
  po [ifi,:,:,:] = p [:,:,:]
  qo [ifi,:,:,:] = q [:,:,:]
  pof[ifi,:,:,:] = pf[:,:,:]
  tpo[ifi] = 1.0*ifi

f.close()
print ("nb points in file: {0}".format(nb_files))

