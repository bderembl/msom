#!/usr/bin/env python

# basic power spectra

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
import fftlib as myfftlib


plt.ion()

def radial_profile(data, center=None):
    y, x = np.indices((data.shape))
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

Lt = 2000
N = 512

Delta = Lt/N
x = np.linspace(0.5*Delta, Lt-0.5*Delta,N)

xx,yy = np.meshgrid(x,x)

kk = 0.02
psi = np.sin(2*np.pi*kk*x)
psi2 = np.sin(2*np.pi*kk*xx)*np.sin(2*np.pi*kk*yy)

#1d
F1 = np.fft.fft(psi)*Delta
F1s = np.fft.fftshift( F1 )
psd1D = np.abs( F1s )**2
kx = np.fft.fftshift(np.fft.fftfreq(N,Delta)) 
dk = kx[1] - kx[0]

#2d
F2 = np.fft.fft2(psi2)*Delta**2
F2s = np.fft.fftshift( F2 )
psd2D = np.abs( F2s )**2
psd2D_r = radial_profile(psd2D)

kr,spec_1D = myfftlib.radial_average(psd2D,Delta)
kr2,spec_1D2 = myfftlib.get_spec_1D(psi2,psi2,Delta)
psd_tmp = myfftlib.get_spec_2D(psi2,psi2,Delta)

# parseval: E = np.sum(psi**2)*Delta = np.sum(F1*F1.conj())*dk
# parseval: E = np.sum(psi2**2)*Delta**2 = np.sum(F2*F2.conj()).real*dk**2 ~= np.sum(spec_1D)*dk*2*np.pi


#kx = np.fft.fftshift(np.fft.fftfreq(N,Delta)) 

# nmax = int(len(psd1D)/np.sqrt(2))
# psd1D = psd1D[:nmax]
# kr = kx[-nmax:]
    

# psi = np.zeros((N,N))
# psi[100:105,:] = 10

# kr = radial_profile(psi)

plt.figure()
#plt.loglog(kx,psd2D/N/Delta,'k-',label=r'PE')
plt.plot(kx,psd1D,'k-',label=r'PE')
#plt.loglog(kr,psd1D,'k-',label=r'PE')
