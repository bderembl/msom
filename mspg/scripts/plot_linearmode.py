#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import interpolate
import matplotlib as mpl
import scipy.io.netcdf as netcdf
from scipy.ndimage.filters import gaussian_filter
import matplotlib

plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


flag_savefig = 1

flag_savenc = 1

# phys consts
L = 5e3

file2 = 'stability_data_4l.npy'

save_dat = np.load(file2,allow_pickle=True)

# with open(file2, "rb") as fp:
#   save_dat = pickle.load(fp)
      
si_d = len(save_dat)
si_y = int(np.sqrt(si_d))
si_x = si_y

# physical grid
xx = np.linspace(0,1,si_x)
yy = np.linspace(0,1,si_y)

xg,yg = np.meshgrid(xx,yy)

x_comp  = np.zeros((si_y*si_x))
y_comp  = np.zeros((si_y*si_x))
o_grid  = np.zeros((si_y*si_x))
k_grid  = np.zeros((si_y*si_x))
l_grid  = np.zeros((si_y*si_x))
x_ind   = np.zeros((si_y*si_x),dtype='int')
y_ind   = np.zeros((si_y*si_x),dtype='int')
Rd_grid = np.zeros((si_y*si_x))

# uncrompress data
nco = 0
for ny in range(0,si_y):
  for nx in range(0,si_x):
    try:
      if save_dat[ny + si_y*nx]['conv'] == 1:
        tmp = save_dat[ny + si_y*nx]['omax']
        if tmp > 0.0:
          x_ind[nco]   = nx
          y_ind[nco]   = ny
          x_comp[nco]  = xx[nx]
          y_comp[nco]  = yy[ny]
          o_grid[nco]  = tmp 
          k_grid[nco]  = save_dat[ny + si_y*nx]['kmax']
          l_grid[nco]  = save_dat[ny + si_y*nx]['lmax']
          Rd_grid[nco] = save_dat[ny + si_y*nx]['Rd'][1] # in km
          
          nco = nco + 1
    except TypeError:
      tmp = 0.0

x_comp  = x_comp[:nco]
y_comp  = y_comp[:nco]
o_grid  = o_grid[:nco]
k_grid  = k_grid[:nco]
l_grid  = l_grid[:nco]
Rd_grid = Rd_grid[:nco]

k_norm = np.sqrt(k_grid**2 + l_grid**2)

# interpolate
o_grid2  = interpolate.griddata((x_comp, y_comp), o_grid,  (xg, yg), method='linear')
k_grid2  = interpolate.griddata((x_comp, y_comp), k_grid,  (xg, yg), method='linear')
l_grid2  = interpolate.griddata((x_comp, y_comp), l_grid,  (xg, yg), method='linear')
k_norm3  = interpolate.griddata((x_comp, y_comp), k_norm,  (xg, yg), method='linear')
Rd_grid2 = interpolate.griddata((x_comp, y_comp), Rd_grid, (xg, yg), method='linear')

# fix nan interpolation
o_grid2_n  = interpolate.griddata((x_comp, y_comp), o_grid,  (xg, yg), method='nearest')
k_grid2_n  = interpolate.griddata((x_comp, y_comp), k_grid,  (xg, yg), method='nearest')
l_grid2_n  = interpolate.griddata((x_comp, y_comp), l_grid,  (xg, yg), method='nearest')
k_norm3_n  = interpolate.griddata((x_comp, y_comp), k_norm,  (xg, yg), method='nearest')
Rd_grid2_n = interpolate.griddata((x_comp, y_comp), Rd_grid, (xg, yg), method='nearest')

o_grid2  = np.where(np.isnan(o_grid2),o_grid2_n,o_grid2)
k_grid2  = np.where(np.isnan(k_grid2),k_grid2_n,k_grid2)
l_grid2  = np.where(np.isnan(l_grid2),l_grid2_n,l_grid2)
k_norm3  = np.where(np.isnan(k_norm3),k_norm3_n,k_norm3)
Rd_grid2 = np.where(np.isnan(Rd_grid2),Rd_grid2_n,Rd_grid2)

lam = 2*np.pi/k_norm3

#o_grid2 = np.where(o_grid2 <= 0, np.NaN, o_grid2)
# in days^-1
o_grid2 = o_grid2*86400.0

# in km^-1
k_grid2 = k_grid2*1e3
l_grid2 = l_grid2*1e3

k_norm2 = np.sqrt(k_grid2**2 + l_grid2**2)

nb_cont = 40

plt.figure()
timescale = 1/o_grid2
maxtimesc = 200.0
timescale2 = np.where(timescale>maxtimesc,maxtimesc,timescale)
timescale2 = np.where(timescale2<1e-1,1e-1,timescale2)
cf = plt.contourf(xx*L,yy*L,o_grid2,nb_cont,cmap=plt.cm.hot_r)
for c in cf.collections:
  c.set_edgecolor("face")

cbar = plt.colorbar(label=r'day$^{-1}$',ticks=[0.1, 0.2, 0.3, 0.4, 0.5])
# cbar = plt.colorbar(label=r'day$^{-1}$',ticks=[0.05, 0.2, 0.5])
# cbar.ax.set_yticklabels([r'$\frac{1}{20}$', r'$\frac{1}{5}$', r'$\frac{1}{2}$'])
plt.contour(xx*L,yy*L,timescale,[-1000,25,100,1e10],colors='k')
plt.contourf(xx*L,yy*L,timescale,[-1000,25,100,1e10],colors='none',hatches=['.', '/', None, '\\\\', '*'],  extend='lower')
plt.text(2500,4100,'1')
plt.text(2500,2050,'2')
plt.text(2500,300,'3')
plt.xlim(0,L)
plt.ylim(0,L)
plt.xlabel('x (km)')
plt.ylabel('y (km)')

if flag_savefig:
  plt.savefig('omax.pdf',bbox_inches='tight')

plt.figure()
cf = plt.contourf(xx*L,yy*L,k_norm3*1e3,nb_cont,cmap=plt.cm.hot_r)
for c in cf.collections:
  c.set_edgecolor("face")

plt.colorbar(label=r'km$^{-1}$',ticks=[0.05, 0.1, 0.15, 0.2])

plt.contour(xx*L,yy*L,timescale,[-1000,25,100,1e10],colors='k')
plt.contourf(xx*L,yy*L,timescale,[-1000,25,100,1e10],colors='none',hatches=['.', '/', None, '\\\\', '*'],  extend='lower')


## define sigma filter
omax = np.nanmax(o_grid2)

omax = 0.25
sigma = o_grid2/omax
sigma = np.where(sigma>1,1,sigma)

sigma = sigma*10 + (1-sigma)*60
sigma_f = gaussian_filter(sigma,5.0)

# plt.contourf(xx*L,yy*L,sigma)
# plt.colorbar(label='km')

plt.xlim(0,L)
plt.ylim(0,L)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
if flag_savefig:
  plt.savefig('k2_omax.pdf',bbox_inches='tight')

plt.figure()
cf = plt.contourf(xx*L,yy*L,2*np.pi/(Rd_grid2*k_norm3),nb_cont,cmap=plt.cm.hot_r)
for c in cf.collections:
  c.set_edgecolor("face")

plt.colorbar(ticks=[2.5, 5,7.5,10,12.5])
plt.xlim(0,L)
plt.ylim(0,L)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
if flag_savefig:
  plt.savefig('k2Rd_omax.pdf',bbox_inches='tight')



fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('click on points')

psi = 1.0*k_norm

plt.scatter(x_comp,y_comp,c=psi, picker=5,cmap=plt.cm.hot_r)


def onpick(event):
  ind = event.ind

  nx = np.take(x_ind, ind)
  ny = np.take(y_ind, ind)
  nxy = ny[0] + si_y*nx[0]
  Rd =      save_dat[nxy]['Rd']
  kt =      save_dat[nxy]['kg']
  lt =      save_dat[nxy]['lg']
  omegai_c =  save_dat[nxy]['omega']
  ktg,ltg = np.meshgrid(kt,lt)
  Rd1 = Rd[1]
  
  plt.figure()
  plt.contourf(ktg*Rd1,ltg*Rd1,(omegai_c[:,:].real*86400),20)
  plt.contourf(ktg*Rd1,-ltg*Rd1,(omegai_c[:,::-1].real*86400),20)
  plt.colorbar(label=r'$\omega (day^{-1}$)')
  plt.xlabel('k x Rd/2pi')
  plt.ylabel('l x Rd/2pi')

  print ('onpick scatter:', ind, np.take(x_comp, ind), np.take(y_comp, ind), np.take(psi, ind), Rd1)
    
fig.canvas.mpl_connect('pick_event', onpick)

plt.show()

# smooth

if flag_savenc:
  fileout = 'sigma.nc'

  f3 = netcdf.netcdf_file(fileout,'w')
  
  f3.createDimension('ypo',si_y)
  f3.createDimension('xpo',si_x)
    
  ypo = f3.createVariable('ypo', 'd', ('ypo',))
  xpo = f3.createVariable('xpo', 'd', ('xpo',))
  
  sigma_o = f3.createVariable('sigma' , 'd', ('ypo','xpo',))

  ypo[:] = np.arange(si_y)
  xpo[:] = np.arange(si_x)
  
  sigma_o [:,:] = sigma_f
  
  f3.close()
