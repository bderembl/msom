#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import interpolate
import matplotlib as mpl
from scipy.ndimage.filters import gaussian_filter
import matplotlib

plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['hatch.linewidth'] = 0.1

flag_savefig = 0
flag_adim_local = 0 #0: adim with l_qg, 1: adim with local deformation radius

# phys consts
L = 100
# qg scales
u_qg = 0.1  # m/s
l_qg = 50e3 # m

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
k_norm2  = interpolate.griddata((x_comp, y_comp), k_norm,  (xg, yg), method='linear')
Rd_grid2 = interpolate.griddata((x_comp, y_comp), Rd_grid, (xg, yg), method='linear')

# fix nan interpolation
o_grid2_n  = interpolate.griddata((x_comp, y_comp), o_grid,  (xg, yg), method='nearest')
k_grid2_n  = interpolate.griddata((x_comp, y_comp), k_grid,  (xg, yg), method='nearest')
l_grid2_n  = interpolate.griddata((x_comp, y_comp), l_grid,  (xg, yg), method='nearest')
k_norm2_n  = interpolate.griddata((x_comp, y_comp), k_norm,  (xg, yg), method='nearest')
Rd_grid2_n = interpolate.griddata((x_comp, y_comp), Rd_grid, (xg, yg), method='nearest')

o_grid2  = np.where(np.isnan(o_grid2),o_grid2_n,o_grid2)
k_grid2  = np.where(np.isnan(k_grid2),k_grid2_n,k_grid2)
l_grid2  = np.where(np.isnan(l_grid2),l_grid2_n,l_grid2)
k_norm2  = np.where(np.isnan(k_norm2),k_norm2_n,k_norm2)
Rd_grid2 = np.where(np.isnan(Rd_grid2),Rd_grid2_n,Rd_grid2)

#o_grid2 = np.where(o_grid2 <= 0, np.NaN, o_grid2)
if flag_adim_local == 0:
  o_grid2 = o_grid2*l_qg/u_qg
  k_norm2 = k_norm2*l_qg/(2*np.pi)
  file_app = 'glob'
else:
  o_grid2 = o_grid2*Rd_grid2/u_qg
  k_norm2 = k_norm2*Rd_grid2/(2*np.pi)
  file_app = 'loc'

nb_cont = 40

# plot omega
plt.figure()
cf = plt.contourf(xx*L,yy*L,o_grid2,nb_cont,cmap=plt.cm.hot_r)
for c in cf.collections:
  c.set_edgecolor("face")

cbar = plt.colorbar(label=r'$\omega$ (t$^{-1}$)',format='%.1f')
plt.contour(xx*L,yy*L,o_grid2,[-1000,0.1,0.5,1e10],colors='k',linewidths=0.5)
plt.contourf(xx*L,yy*L,o_grid2,[-1000,0.1,0.5,1e10],colors='none',hatches=['x', '/', None, '\\\\', '*'],  extend='lower')
plt.text(70,70,'1')
plt.text(50,40,'2')
plt.text(50,6,'3')
plt.xlabel('x')
plt.ylabel('y')
if flag_savefig:
  plt.savefig('omax_' + file_app + '.pdf',bbox_inches='tight')

# plot k
plt.figure()
cf = plt.contourf(xx*L,yy*L,k_norm2,nb_cont,cmap=plt.cm.hot_r)
for c in cf.collections:
  c.set_edgecolor("face")

plt.colorbar(label=r'$k$ (cycle/l)',format='%.1f')

plt.contour(xx*L,yy*L,o_grid2,[-1000,0.1,0.5,1e10],colors='k',linewidths=0.5)
plt.contourf(xx*L,yy*L,o_grid2,[-1000,0.1,0.5,1e10],colors='none',hatches=['x', '/', None, '\\\\', '*'],  extend='lower')
plt.xlabel('x')
plt.ylabel('y')
if flag_savefig:
  plt.savefig('k2_omax_' + file_app + '.pdf',bbox_inches='tight')

# clickable plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('click on points')

psi = 1.0*k_norm
plt.scatter(x_comp*si_x,y_comp*si_x,c=psi, picker=5,cmap=plt.cm.hot_r)


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
  plt.contourf(ktg*Rd1,ltg*l_qg/(2*np.pi),(omegai_c[:,:].real*l_qg/u_qg),20)
  plt.contourf(ktg*Rd1,-ltg*l_qg/(2*np.pi),(omegai_c[:,::-1].real*l_qg/u_qg),20)
  plt.colorbar(label=r'$\omega (day^{-1}$)')
  plt.xlabel('k (cycle/l)')
  plt.ylabel('l (cycle/l)')

  print ('onpick scatter:', ind, np.take(x_comp, ind), np.take(y_comp, ind), np.take(psi, ind), Rd1)
    
fig.canvas.mpl_connect('pick_event', onpick)

plt.show()
