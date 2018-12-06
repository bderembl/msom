import numpy as np
import matplotlib.pyplot as plt
import glob,os,re

plt.ion()

dir0 = "../outdir/"

file_b = 'b*'
file_u = 'u*'

allfilesu = sorted(glob.glob(dir0 + file_u));
allfilesb = sorted(glob.glob(dir0 + file_b));
nb_files  = len(allfilesu);

# scales
# L (m)
L = 5000e3

#H (M)
H = 5000

# beta (1/(ms))
beta = 2.0e-11

# N2 (1/s**2)
N2 = 1e-6

Ts = beta*L**3/(N2*H**2)
Us = N2*H**2/(beta*L**2)
Ws = N2*H**3/(beta*L**3)
Bs = N2*H
Thetas = Bs/10/2e-4 # 1/g alpha
Kv = N2*H**4/(beta*L**3)
Kh = N2*H**2/(beta*L)

b = np.fromfile(allfilesb[0],'f4')
N = int(b[0])
N1 = N + 1
nl2 = int(len(b)/N1**2)
nl = nl2 - 2


Delta = 1/N
x = np.linspace(0.5*Delta, 1-0.5*Delta,N)
y = np.linspace(0.5*Delta, 1-0.5*Delta,N)

xc, yc = np.meshgrid(x,y)

ifi = nb_files - 1

#temporary
#z = np.linspace(-1,0,nl)
z = np.linspace(0,-1,nl)

for ifi in range(0,nb_files):
#for ifi in range(nb_files-1,nb_files):

  b = np.fromfile(allfilesb[ifi],'f4').reshape(nl2,N1,N1).transpose(0,2,1)
  uv = np.fromfile(allfilesu[ifi],'f4').reshape(2*nl2,N1,N1).transpose(0,2,1)

  b = b[1:-1,1:,1:]
  u = uv[2:-2:2,1:,1:]
  v = uv[3:-2:2,1:,1:]


  # choose vertical level (0=top, nl-1=bottom)
  l = 0
  
  #plot velocity vector every nsk pts
  nsk = 2
  
  vmin = np.min(b)
  vmax = np.max(b)
  
  plt.figure(1)
  iy = int(N/2)-1
  plt.clf()
  plt.subplot(211)
  iy = int(N/2)-1
  plt.contourf(x,z,b[:,iy,:],15,vmin=vmin,vmax=vmax)
  plt.contour(x,z,b[:,iy,:],15,vmin=vmin,vmax=vmax,colors='w',linewidths=0.5)
  #plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("z")
  
  plt.subplot(212)
  ix = int(N/2)-1
  plt.contourf(y,z,b[:,:,ix],15,vmin=vmin,vmax=vmax)
  plt.contour(y,z,b[:,:,ix],15,vmin=vmin,vmax=vmax,colors='w',linewidths=0.5)
  #plt.colorbar()
  plt.xlabel("y")
  plt.ylabel("z")
  plt.draw()
  
  # plt.pcolormesh(x[0,:,:],y[0,:,:],b[l,:,:],vmin=vmin,vmax=vmax)
  # plt.colorbar()
  # plt.contour(x[0,:,:],y[0,:,:],topo,colors='w',linewidths=0.5)
  # plt.xlabel('x')
  # plt.ylabel('y')
  
  
  plt.figure(2)
  plt.clf()
  plt.contourf(xc,yc,b[l,:,:],25, cmap=plt.cm.RdYlBu_r)
#  plt.pcolormesh(x,y,w[:,:,l])
  plt.colorbar()
#  plt.contour(x,y,topo[:,:,0],colors='w',linewidths=0.5)
#  plt.contour(x,y,psibt,10,colors='r',linewidths=0.5)
  plt.quiver(xc[::nsk,::nsk],yc[::nsk,::nsk],u[l,::nsk,::nsk],v[l,::nsk,::nsk])
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Buoyancy (color), velocity (arrows), topography (contours), and BT stream function (red)')
  plt.draw()
#  plt.savefig('xyplot.png')
  
  

  plt.figure(3)
  plt.clf()
#  iy = np.argmin(np.abs(y[0,:]-1.1))
  psi = v[:,iy,:]
  vmax = np.max(np.abs(psi))
  plt.contour(x,z,b[:,iy,:],colors="k",linewidths=1.0)
  plt.contourf(x,z,psi,30,cmap=plt.cm.bwr, vmax=vmax, vmin=-vmax)
  # plt.plot(x[:,iy],z[:,iy,0],"k",linewidth=1.5)
  # plt.fill_between(x[:,iy], z[:,iy,0],np.min(z[:,iy,0]) + 0*z[:,iy,0],facecolor="0.7")
  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("z")
  plt.title('V vel (color), buoyancy (contours)')
  plt.draw()
#  plt.savefig('xzplot.png')


#   plt.figure(1)
#   plt.clf()
# #  plt.imshow(u[:,:,5])
#   plt.plot(b[15,15,:])
# #  plt.colorbar()
#   plt.draw()
  input("maxu = {0:.2f}, maxv = {1:.2f}".format(Us*np.max(u[l,:,:]),Us*np.max(v[l,:,:])))
  
