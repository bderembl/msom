#!/usr/bin/env python

import scipy.io.netcdf as netcdf
import numpy as np
import matplotlib.pyplot as plt
import def_radius
import sympy as sp
import scipy.linalg as la
import scipy.optimize
import os.path
import matplotlib

plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


filepsi = '../psipg.bas'
filegp  = '../gppg.bas'
fileh   = '../dh.bin'

file2 = 'stability_data_4.npy'

flag_eigmeth = 1 # 0: via characteristic polynomial (slower?), 1: eigenvalue solver
flag_test2lay = 0
flag_gradqbar = 0 # format of the large scale advection of PV 0: L grad psi, 1: grad L psi
flag_plot = 0
flag_print = 1
flag_klmap = 1
flag_savedat = 1

if flag_savedat:
  if os.path.isfile(file2):
    save_dat = np.load(file2)
  else:
    save_dat = []
else:
  save_dat = []
  

# dimensions
b = np.fromfile(filepsi,'f4')
N = int(b[0])
N1 = N + 1
nl = int(len(b)/N1**2)

# PG scales
L = 5000e3  # m
H = 5000    # m
beta = 2.0e-11 # 1/m/s
N2 = 1e-6  #  (1/s**2)

Ts = beta*L**3/(N2*H**2)
Us = N2*H**2/(beta*L**2)
Ws = N2*H**3/(beta*L**3)
Bs = N2*H
Thetas = Bs/10/2e-4 # 1/g alpha
Kv = N2*H**4/(beta*L**3)
Kh = N2*H**2/(beta*L)
ys = 0.3

# grid
Delta = 1/N
Deltad = L*Delta
x = np.linspace(0.5*Delta, 1-0.5*Delta,N)
y = ys + np.linspace(0.5*Delta, 1-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

f0 = yc*L*beta


# read data
gp = np.fromfile(filegp,'f4').reshape(nl,N1,N1).transpose(0,2,1)
po = np.fromfile(filepsi,'f4').reshape(nl,N1,N1).transpose(0,2,1)
th = np.fromfile(fileh,'f4')

gp = gp[:,1:,1:]
psib = po[:,1:,1:]

##### TESTING PART
# make a 2 layer problem
if flag_test2lay:
  nl = 2
  po = po[:nl,:,:]
  psib = psib[:nl,:,:]
  
  th = th[:nl]
  gp = gp[:nl,:,:]
  
  
  
  # #idealized profile
  f0 = 1e-4 + 0*f0
  u1 =  0.0001
  u2 = -0.0004
  u3 =  0.01
  du = u2-u1

  psib[0,:,:] = -u1*yc*L*f0
  psib[1,:,:] = -u2*yc*L*f0
  # #psib[2,:,:] = -u3*yg
  # #psib[3:,:,:] = 0.001*yg
  
  th = 0*th + 200.0
  gp = 0*gp + 1.0e-4
##### end testing part


aux,dpsibdy,dpsibdx = np.gradient(psib)
dpsibdy = dpsibdy/Deltad
dpsibdx = dpsibdx/Deltad

for nz in range(0,nl):
  dpsibdy[nz,:,:] = dpsibdy[nz,:,:]/f0
  dpsibdx[nz,:,:] = dpsibdx[nz,:,:]/f0
  
dqbdy = np.zeros((nl,N,N))
dqbdx = np.zeros((nl,N,N))
for nx in range(0,N):
  for ny in range(0,N):
    mata = def_radius.construct_mat(th,gp[:,ny,nx],f0[ny,nx])
    dqbdy[:,ny,nx] = -np.dot(mata,dpsibdy[:,ny,nx])
    dqbdx[:,ny,nx] = -np.dot(mata,dpsibdx[:,ny,nx])
      
# add beta effect
dqbdy = dqbdy + beta

if len(save_dat) == 0:
  save_dat = [None]*N*N
  
# # specific location
# figure middle domain: 320,320
nx1 = 32
ny1 = 32
nx2 = nx1+1
ny2 = ny1+1
nsk = 1

nx1 = 0
ny1 = 0
nx2 = N
ny2 = N
nsk = 1

k_grid = np.zeros((N,N))
l_grid = np.zeros((N,N))
o_grid = np.zeros((N,N))

ntot = len(range(nx1,nx2,nsk))*len(range(ny1,ny2,nsk))
if nx2-nx1 > 1 or ny2-ny1 > 1:
  print ('computing', ntot, 'points')
nco = 0
nco1 = 0
for nx in range(nx1,nx2,nsk):
  for ny in range(ny1,ny2,nsk):
    nco1 = nco1 + 1

    # periodic save
    if nco%20 == 0 and flag_savedat:
      np.save(file2, save_dat)

    try:
      save_dat[ny + N*nx]['conv'] == 1
      print('oldconv: ', save_dat[ny + N*nx]['conv'])
      if save_dat[ny + N*nx]['conv'] > 0:
        doit = 0
      else:
        doit = 2
    except TypeError:
      doit = 1
    if doit > 0: 
      nco = nco + 1
      if nx2-nx1 > 1 or ny2-ny1 > 1:
        print (nco,'; ',nco1, '/', ntot, '  ', ny + N*nx)
      # compute deformation radius
      Rd = def_radius.cal_rad(th,gp[:,ny,nx],f0[ny,nx])
      #Rd = def_radius.cal_rad(th,gp[:,ny,nx],fmid)
      Rd1 = Rd[1]
      
      mata = def_radius.construct_mat(th,gp[:,ny,nx],f0[ny,nx])
      #mata = def_radius.construct_mat(th,gp[:,ny,nx],fmid)

      l2m,m2l = def_radius.cal_transfo(th,gp[:,ny,nx],f0[ny,nx])
      
      # compute matrices
      k,l,o,lam = sp.symbols('k l o lam')
      k2 = k**2 + l**2
      idz = sp.eye(nl)
      idznp = np.eye(nl)
      
      p2q = -mata-k2*idz
      p2qnp = -mata-k2*idznp
            
      
      diag1 = k*dqbdy[:,ny,nx] - l*dqbdx[:,ny,nx]
      diag2 = l*dpsibdx[:,ny,nx] - k*dpsibdy[:,ny,nx]
      
      diag1 = np.diagflat(diag1)
      diag2 = np.diagflat(diag2)
      
      mat1 = diag1 + np.dot(diag2,p2qnp)
      
      if flag_eigmeth == 0:
        ip2q = p2q.inv()
        mat2 = ip2q*sp.Matrix((mat1))
        print ('computing characteristic polynomial')
        poly = mat2.charpoly(lam)
        print ('.. done')
      
      
      Rd2 = Rd1**2
      # si_k = 75
      # si_l = 37

      # figure middle domain: 150,75
      si_k = 50
      si_l = 25

      # figure middle of the domain
      #modes = 2*np.array([[48,8],[38,8],[43,2]],np.int32)

      modes = 2*np.array([[48,8]],np.int32)
      nmodes,naux = modes.shape
      
      # size of the k-l window
      if doit == 1:
      # figure middle domain: 3.5
        window_s = 3.5
#        window_s = 30
      else:
        window_s = save_dat[ny + N*nx]['window_s']
        
      k0min = -window_s/(Rd1)
      k0max = window_s/(Rd1)
      # get l<0 by symmetry
      l0min = 1e-3/(Rd1)
      l0max = window_s/(Rd1)
      kt = np.linspace(k0min,k0max,si_k)
      lt = np.linspace(l0min,l0max,si_l)
      
      ktg,ltg = np.meshgrid(kt,lt)
      
      if flag_test2lay:
        omegai_th = np.imag(-0.5*du*kt*(1./((Rd1**2*kt**2)*(Rd1**2*kt**2+1.)))*np.sqrt(beta**2*Rd1**4/du**2-Rd1**4*kt**4*(1.-Rd1**4*kt**4)+0j))
        # simplified version without beta
        omegai_thwob = np.imag(-0.5*du*kt*np.sqrt(-(1.-Rd2*kt**2)/(1.+Rd2*kt**2)+0j))
        #omegai_thwob = -0.5*du*kt*np.sqrt((1.-Rd2*kt**2)/(1.+Rd2*kt**2))
        omegar_th = (u2+u1)*0.5*kt
      
      omegai_c = np.zeros((si_l,si_k,nl),dtype=complex)
      eivec_c  = np.zeros((si_l,si_k,nl,nl),dtype=complex)
      
      
      if flag_klmap:
        # ke = 0.1/(Rd[1]*1e3)
        # #le = 1/(Rd[1]*1e3)
        # le = 0.0 #0.1/(Rd[1]*1e3)
        for il in range(0,si_l):
          for ik in range(0,si_k):
            ke = kt[ik]
            le = lt[il]
            if flag_eigmeth == 0:
              poly2 = poly.subs([(k,ke),(l,le)])
              omega = sp.roots(poly2,lam).items()
              nv = len(omega)
              for n in range(0,nv):
                omegai_c[il,ik,n] = omega[n][0]
          
              omegai_c[il,ik,:nv] = np.sort(np.imag(omegai_c[il,ik,:nv]))[::-1]
          
            elif flag_eigmeth == 1:       # # other method (generalized eigenvalue)
          
              mat1_s = np.float32(np.array(sp.Matrix((mat1)).subs([(k,ke),(l,le)])))
              mat3_s = np.float32(np.array(p2q.subs([(k,ke),(l,le)])))
              # not sure which one is faster for eigenvalues only (seems about the same time)
              eival,eivec = la.eig(mat1_s,mat3_s)
              #eival = la.eigvals(mat1_s,mat3_s)
            
              #omegai_c[il,ik,:] = eival
              idx = np.argsort(np.imag(eival))[::-1]
              omegai_c[il,ik,:] = np.imag(eival[idx])
              eivec_c[il,ik,:,:] = np.real(eivec[:,idx])
              
        (ilmax,ikmax) = np.unravel_index(omegai_c[:,:,0].argmax(), omegai_c[:,:,0].shape)
        #ikmax = np.argmax(omegai_c[:,0])

        
        kmax = kt[ikmax]
        lmax = lt[ilmax]
        omax = np.real(omegai_c[ilmax,ikmax,0])

        eivecmax = eivec_c[ilmax,ikmax,:,0].squeeze()
        
        fw = np.int(si_l/4)
        
        om_bdy0 = np.sum(omegai_c[-1,:,0]) + np.sum(omegai_c[:,0,0]) + np.sum(omegai_c[:,-1,0])
        om_bdy1 = np.sum(omegai_c[-1-fw,:,0]) + np.sum(omegai_c[:,+fw,0]) + np.sum(omegai_c[:,+fw,0])

        myconv = 0
        
        if ikmax > 0 and ikmax < si_k-1:
          k_grid[ny,nx] = kmax
          l_grid[ny,nx] = lmax
          o_grid[ny,nx] = omax
          myconv = 1
          
          if np.abs(om_bdy0) > 1e-7:
            myconv = -1 # enlarge window
            window_s = 1.2*window_s
            print('enlarging window',np.abs(om_bdy0))
          elif np.abs(om_bdy1) < 1e-11:
            myconv = -1 # reduce window
            window_s = 0.8*window_s
            print('reducing window',np.abs(om_bdy0))

          if flag_print:
            print ('Deformation Radius:', Rd1*1e-3)
            print ('kmax x Rd1 = ',kmax*Rd1)
            print ('lmax x Rd1 = ',lmax*Rd1)
            print ('time scale = ',np.real(1/omax/86400), ' days')
            print ('projection on vert modes:', np.dot(l2m[:,:],eivecmax.real))
            
        else:
          if np.sum(omegai_c[:,:,0]) == 0: # no instability
            myconv = 2
          else:
            myconv = -1 # enlarge window
            window_s = 1.5*window_s
            
          k_grid[ny,nx] = np.NaN
          l_grid[ny,nx] = np.NaN
          o_grid[ny,nx] = np.NaN
          if flag_print:
            print ('no reliable unstable mode found -> increase the k-l window')

            
        save_dat[ny + N*nx] = {'conv':myconv, 
                                  'Rd':Rd,
                                  'kmax':kmax,
                                  'lmax':lmax,
                                  'omax': omax,
                                  'kg':kt,
                                  'lg':lt,
                                  'omega':omegai_c[:,:,0].squeeze(),
                                  'window_s': window_s}

      else: #if flag_klmap
        def find_omegamax(xk):
          ke = xk[0]
          le = xk[1]
        
          mat1_s = np.float32(np.array(sp.Matrix((mat1)).subs([(k,ke),(l,le)])))
          mat3_s = np.float32(np.array(p2q.subs([(k,ke),(l,le)])))
          eival = la.eigvals(mat1_s,mat3_s)
            
          omegai_c = eival
          omegai_c = np.sort(np.imag(omegai_c))[::-1]
    #      print ke,le,omegai_c[0]
          # return - because looking for min
          return -omegai_c[0]
        
        theta = np.arctan2(dpsibdx[0,ny,nx],-dpsibdy[0,ny,nx])
        ke = np.cos(theta)/(Rd[1]*1e3)
        le = np.sin(theta)/(Rd[1]*1e3)
        
        [kmax,lmax] = scipy.optimize.fmin_powell(lambda x: find_omegamax(x),[ke,le],ftol=1e-1,disp=False)
        omax = -find_omegamax([kmax,lmax])
        k_grid[ny,nx] = kmax
        l_grid[ny,nx] = lmax
        o_grid[ny,nx] = omax
          
  
      # 2d plot
      if flag_plot and flag_klmap:
        plt.figure()
        plt.contourf(ktg*Rd1,ltg*Rd1,(omegai_c[:,:,0].real*86400),20)
        plt.contourf(ktg*Rd1,-ltg*Rd1,(omegai_c[:,::-1,0].real*86400),20)
        cb = plt.colorbar(label=r'$\omega$ (day$^{-1}$)')
        cb.formatter.set_powerlimits((0, 0))
        cb.update_ticks()
        
        plt.xlabel('k x Rd')
        plt.ylabel('l x Rd')

        modes[0,0] = ikmax
        modes[0,1] = ilmax
        for i in range(0,nmodes):
          ikmax2 = modes[i,0]
          ilmax2 = modes[i,1]
          plt.plot(kt[ikmax2]*Rd1,lt[ilmax2]*Rd1,'ro',markersize=10)
          plt.text(kt[ikmax2]*Rd1,lt[ilmax2]*Rd1,str(i+1),ha='center',va='center')

          
          kmax2 = kt[ikmax2]
          lmax2 = lt[ilmax2]
          omax2 = np.real(omegai_c[ilmax2,ikmax2,0])

          eivecmax2 = eivec_c[ilmax2,ikmax2,:,0].squeeze()
          if flag_print:
            print ('mode:', i+1)
            print ('kmax x Rd1 = ',kmax2*Rd1)
            print ('lmax x Rd1 = ',lmax2*Rd1)
            print ('time scale = ',np.real(1/omax2/86400), ' days')
            print ('projection on vert modes:', np.dot(l2m[:,:],eivecmax2.real))
            
        #plt.savefig('klo_map.png',bbox_inches='tight')
        #plt.savefig('klo_map_340_340.pdf',bbox_inches='tight')
        #plt.savefig('klo_map_340_340_30l_1.pdf',bbox_inches='tight')
        
      # 1d plot
      if flag_test2lay:
        plt.figure()
        plt.plot(kt*Rd1,omegai_c[0,:,0]*86400,'k')
        plt.plot(kt*Rd1,omegai_th*86400,'r--')
        plt.plot(kt*Rd1,omegai_thwob*86400,'b--')
        plt.xlim([0.0,kt[-1]*Rd1])
        
        plt.xlabel('k x Rd')
        plt.ylabel(r'$\omega (day^{-1})$')
  

# if nx2-nx1 > 1 and ny2 - ny1>1:
#   plt.figure()
#   timescale = 1/o_grid[::nsk,::nsk]/86400.0
#   maxtimesc = 200.0
#   timescale2 = np.where(timescale>maxtimesc,maxtimesc,timescale)
#   plt.contourf(xx[::nsk]*1e-3,yy[::nsk]*1e-3,timescale2,100,cmap=plt.cm.hot)
#   plt.colorbar()
#   #plt.savefig('timescale_instab.png',bbox_inches='tight')

if flag_savedat:
  np.save(file2, save_dat)
