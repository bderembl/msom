#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
sys.path.append('../')
import pypg as bas

plt.ion()

N = 64
nl = 4
N1 = N + 1
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X,Y = np.meshgrid(x,y)
ne = N*N*nl + 2*N*N1*nl

bas.init_grid(N)
bas.pyinit_const(nl)
bas.set_vars()
bas.pyinit_last()

var0 = 0*np.random.rand(ne)
dvar = np.zeros(ne)

bas.pystep(var0,dvar)

b0 = dvar[:N*N*nl]
u0 = dvar[N*N*nl:N*N*nl+N1*N*nl]
v0 = dvar[N*N*nl+N1*N*nl:]

b1 = b0.reshape((N,N,nl))
u1 = u0.reshape((N1,N,nl))
v1 = v0.reshape((N,N1,nl))

#bas.pytrash_vars()
#bas.run()

#plt.contourf(X,Y,db0[-1,:,:].T); plt.colorbar()
