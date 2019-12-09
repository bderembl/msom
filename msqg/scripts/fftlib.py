#!/usr/bin/env python

import numpy as np


def radial_average(spec_2D,Delta):
  ''' Compute the azimuthal avearge of the 2D spectrum '''

  KK,kr = get_wavenumber(spec_2D,Delta)
  spec_1D = np.zeros(len(kr))
  for i in range(kr.size):
    kfilt =  (KK>=kr[i] - kr[0]) & (KK<=kr[i])
    Nbin = kfilt.sum()
    spec_1D[i] = (spec_2D[kfilt].sum())*kr[i]/Nbin
  return kr,spec_1D


def get_len_wavenumber(N,Delta):
  kx = np.fft.fftshift(np.fft.fftfreq(N,Delta))
  dk = np.abs(kx[2]-kx[1])
  return (int(kx.max()/dk)-1)


def get_wavenumber(spec_2D,Delta):
  ''' Compute wavenumber and radial wavenumber '''
  N,N2 = np.shape(spec_2D) # assume N = N2
  kx = np.fft.fftshift(np.fft.fftfreq(N,Delta)) # two sided  
  k, l = np.meshgrid(kx,kx)
  KK = np.sqrt(k**2 + l**2)  
  dk = np.abs(kx[2]-kx[1])
  kr = dk*np.arange(1,int(k.max()/dk))
  return KK,kr


def get_spec_2D(psi1,psi2,Delta):
  ''' Compute the 2D power spectrum of the data 

   normalization such that parseval is ok: 
   E = np.sum(psi**2)*Delta**2 
     = np.sum(psi_hat*psi_hat.conj()).real*dk**2 
     ~ np.sum(spec_1D)*dk*2*np.pi '''
  
  psi1_hat = np.fft.fft2(psi1)*Delta**2
  psi2_hat = np.fft.fft2(psi2)*Delta**2
  spec_2D = (psi1_hat*psi2_hat.conj()).real
  spec_2D = np.fft.fftshift(spec_2D)
  return spec_2D


def get_spec_1D(psi1,psi2,Delta):
  ''' Compute the 2D power spectrum of the data '''

  spec_2D = get_spec_2D(psi1,psi2,Delta)
  kr,spec_1D = radial_average(spec_2D,Delta)
  return kr,spec_1D


def get_flux(psi1,psi2,Delta):
  ''' Compute flux'''
  spec_2D = get_spec_2D(psi1,psi2,Delta)

  # kr,spec_1D = radial_average(spec_2D,Delta)
  # flux = -np.cumsum(spec_1D)*dk # integrate from low wavenumbers
  # flux = np.cumsum(spec_1D[::-1])[::-1]*dk # integrate from high wavenumbers

  KK,kr = get_wavenumber(spec_2D,Delta)
  dk = kr[1] - kr[0]

  flux = np.zeros(len(kr))
  for i in range(kr.size):
    kfilt =  (kr[i] <= KK ) 
    flux[i] = (spec_2D[kfilt]).sum()*dk*dk
  return kr,flux
  
