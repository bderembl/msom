#!/usr/bin/env python

import numpy as np

N = 128
nl = 4
nptr = 2

Delta = 1/N
x  = np.linspace(0.5*Delta, 1-0.5*Delta,N)

# ptracers are (z,y,x)
# frist ptracer
ptr0 = np.zeros((nl,N,N))
ptr0[0,:,:] = np.cos(2*np.pi*x)
ptr0[1,:,:] = 2*x

# second ptracer
ptr1 = np.zeros((nl,N,N))
ptr1[0,:,:] = np.sin(2*np.pi*x)
ptr1[1,:,:] = 10*x

# combine and output in .bas format: (z,x,y)
ptr = np.zeros((nl*nptr,N+1,N+1))
ptr[:,0,0] = N
ptr[0::nptr,1:,1:] = ptr0
ptr[1::nptr,1:,1:] = ptr1
ptr = np.transpose(ptr,(0,2,1))
ptr.astype('f4').tofile('ptr0.bas')
