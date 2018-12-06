import numpy as np
import matplotlib.pyplot as plt
import glob
import array

plt.ion()

dir0 = "../outdir/"

filep = 'po*'


lref = 50e3
uref = 0.1

allfilesp = sorted(glob.glob(dir0 + filep));
nb_files  = len(allfilesp);

b = np.fromfile(allfilesp[0],'f4')
N = int(b[0])
N1 = N + 1
nl = int(len(b)/N1**2)


y = np.arange(N)
x = np.arange(N)

xc,yc = np.meshgrid(x,y)

psipg = np.fromfile(dir0 + 'psipg.bas0512','f4').reshape(nl,N1,N1).transpose(0,2,1)
psipg = psipg[:,1:,1:]

gppg = np.fromfile(dir0 + 'gppg.bas0512','f4').reshape(nl,N1,N1).transpose(0,2,1)
gppg = gppg[:,1:,1:]


# read constants
iBu = np.fromfile(dir0 + 'iBu.bas','f4').reshape(nl,N1,N1).transpose(0,2,1)
iBu = iBu[:,1:,1:]

Rd = np.sqrt(-1/iBu)*lref

plt.figure()
CS = plt.contour(xc,yc,Rd[1,:,:]*1e-3,[1,2,5,10,20,30,40,50,60,70,80,90],colors='b')
plt.clabel(CS, inline=1, fontsize=10)

# choose vertical level (0= top, nl-1=bottom)
l = 0

#ifi0 = nb_files - 1
ifi0 = 0

for ifi in range(ifi0,nb_files):

  p  = np.fromfile(allfilesp[ifi],'f4').reshape(nl,N1,N1).transpose(0,2,1)
  p = p[:,1:,1:]

  plt.figure(2)
  plt.clf()

  plt.contour(xc,yc,psipg[l,:,:] + p[l,:,:],colors='k',linewidths=1)

  plt.draw()
  
  nextit = input("press enter")
