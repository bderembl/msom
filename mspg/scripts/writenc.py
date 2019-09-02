#!/usr/bin/env python

import numpy as np
import sys,glob,os,re
import scipy.io.netcdf as netcdf

dir0 = "../outdir_"

if len(sys.argv) > 1:
  dir0 = dir0 + str(format(sys.argv[1])).zfill(4) + '/'

exec(open(dir0 + "params.in").read())

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

tpo = f.createVariable('t', 'f', ('t',))
zpo = f.createVariable('z', 'f', ('z',))
ypo = f.createVariable('y', 'f', ('y',))
xpo = f.createVariable('x', 'f', ('x',))

bo  = f.createVariable('b' , 'f', ('t','z','y','x',))
uo  = f.createVariable('u' , 'f', ('t','z','y','x',))
vo  = f.createVariable('v' , 'f', ('t','z','y','x',))

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
  uo[ifi,:,:,:] = 1.*u[:,:,:]
  vo[ifi,:,:,:] = 1.*v[:,:,:]

f.close()
print ("nb points in file: {0}".format(nb_files))

