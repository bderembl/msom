#!/usr/bin/env python

import numpy as np
import glob,os,re
import scipy.io.netcdf as netcdf

dir0 = "../outdir/"

fileb = 'b*'
fileu = 'u*'

allfilesb = sorted(glob.glob(dir0 + fileb));
allfilesu = sorted(glob.glob(dir0 + fileu));
nb_files  = len(allfilesb);

# dimensions
b = np.fromfile(allfilesb[0],'f4')
N = int(b[0])
N1 = N + 1
nl2 = int(len(b)/N1**2)
nl = nl2 - 2

# create netcdf file
fileout = dir0 + 'vars.nc'
f = netcdf.netcdf_file(fileout,'w')

f.createDimension('t',None)
f.createDimension('z',nl)
f.createDimension('y',N)
f.createDimension('x',N)

tpo = f.createVariable('t', 'd', ('t',))
zpo = f.createVariable('z', 'd', ('z',))
ypo = f.createVariable('y', 'd', ('y',))
xpo = f.createVariable('x', 'd', ('x',))

bo  = f.createVariable('b' , 'd', ('t','z','y','x',))

zpo[:] = np.arange(nl)
ypo[:] = np.arange(N)
xpo[:] = np.arange(N)
tpo[:] = np.arange(nb_files)

ifi = nb_files - 1

for ifi in range(0,nb_files):

  b = np.fromfile(allfilesb[ifi],'f4').reshape(nl2,N1,N1).transpose(0,2,1)
  uv = np.fromfile(allfilesu[ifi],'f4').reshape(2*nl2,N1,N1).transpose(0,2,1)

  b = b[1:-1,1:,1:]
  u = uv[2:-2:2,1:,1:]
  v = uv[3:-2:2,1:,1:]
  
  bo[ifi,:,:,:] = 1.*b[:,:,:]

f.close()
print ("nb points in file: {0}".format(nb_files))

