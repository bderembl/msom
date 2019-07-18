#!/usr/bin/env python

# basic power spectra

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack


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


F1 = np.fft.fft(psi)
F2 = np.fft.fftshift( F1 )
psd2D = np.abs( F2 )**2
kx = np.fft.fftshift(np.fft.fftfreq(N,Delta)) 


# F1 = np.fft.fft2(psi2)
# F2 = np.fft.fftshift( F1 )
# psd2D = np.abs( F2 )**2
# psd1D = radial_profile(psd2D)

#kx = np.fft.fftshift(np.fft.fftfreq(N,Delta)) 

# nmax = int(len(psd1D)/np.sqrt(2))
# psd1D = psd1D[:nmax]
# kr = kx[-nmax:]
    

# psi = np.zeros((N,N))
# psi[100:105,:] = 10

# kr = radial_profile(psi)

plt.figure()
plt.loglog(kx,psd2D,'k-',label=r'PE')
#plt.loglog(kr,psd1D,'k-',label=r'PE')
