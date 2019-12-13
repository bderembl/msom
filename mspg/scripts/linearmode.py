#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import def_radius
import scipy.linalg as la
import scipy.optimize
import os.path
import matplotlib
from skimage.feature import peak_local_max
from scipy.ndimage.filters import gaussian_filter

from scipy import sparse

plt.ion()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

nl = 30

filefr = 'frpg_' + str(nl) +'l.bas'
filepsi = 'psipg_' + str(nl) +'l.bas'
fileh = 'dh_' + str(nl) +'l.bin'

file2 = 'stability_data_' + str(nl) + 'l.npy'

flag_test2lay = 0
flag_gradqbar = 0 # format of the large scale advection of PV 0: L grad psi, 1: grad L psi
flag_plot = 1
flag_print = 1
flag_klmap = 1
flag_savedat = 0
flag_dim = 1 #
flag_localRd = 0 # rescale by local Rd (flag=1) or global Rd (Flag=0)

if flag_savedat:
  if os.path.isfile(file2):
    save_dat = np.load(file2,allow_pickle=True)
  else:
    save_dat = []
else:
  save_dat = []
  

# dimensions
b = np.fromfile(filepsi,'f4')
N = int(b[0])
N1 = N + 1
#nl = int(len(b)/N1**2)

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

# qg scales
u_qg = 0.1  # m/s
l_qg = 50e3 # m

# friction
Re  = 50
Re4 = 1000
Ek  = 0.05

# dimensional
nu  = u_qg*l_qg/Re
nu4 = u_qg*l_qg**3/Re4
bf  = Ek*u_qg/l_qg

# grid
Delta = 1/N
Deltad = L*Delta
x = np.linspace(0.5*Delta, 1-0.5*Delta,N)
y = ys + np.linspace(0.5*Delta, 1-0.5*Delta,N)
xc, yc = np.meshgrid(x,y)

f0 = yc*L*beta

# read data
fr = np.fromfile(filefr,'f4').reshape(nl,N1,N1).transpose(0,2,1)
po = np.fromfile(filepsi,'f4').reshape(nl,N1,N1).transpose(0,2,1)
dh = np.fromfile(fileh,'f4')

dh = dh*H
dhi = 0.5*(dh[:-1] + dh[1:])
hh = np.cumsum(dh)

fr = fr[:,1:,1:]
psib = po[:,1:,1:]*u_qg*l_qg
#psib = po[:,1:,1:]*u_qg*u_qg/f0

gp = 1/fr[:-1,:,:]**2*u_qg**2/H**2*dhi.reshape(nl-1,1,1)

##### TESTING PART
# make a 2 layer problem
if flag_test2lay:
  nl = 2
  po = po[:nl,:,:]
  psib = psib[:nl,:,:]
  
  dh = dh[:nl]
  gp = gp[:nl,:,:]
  
  
  
  # #idealized profile
  f0 = 1e-4 + 0*f0
  u1 =  0.0001
  u2 = -0.0004
  u3 =  0.01
  du = u2-u1

  #psib[0,:,:] = -u1*yc*L*f0
  #psib[1,:,:] = -u2*yc*L*f0
  psib[0,:,:] = -u1*yc*L
  psib[1,:,:] = -u2*yc*L
  # #psib[2,:,:] = -u3*yg
  # #psib[3:,:,:] = 0.001*yg
  
  dh = 0*dh + 200.0
  gp = 0*gp + 1.0e-4
##### end testing part


aux,dpsibdy,dpsibdx = np.gradient(psib,Deltad)

  
dqbdy = np.zeros((nl,N,N))
dqbdx = np.zeros((nl,N,N))
for nx in range(0,N):
  for ny in range(0,N):
    mata = def_radius.construct_mat(dh,gp[:,ny,nx],f0[ny,nx])
    dqbdy[:,ny,nx] = -np.dot(mata,dpsibdy[:,ny,nx])
    dqbdx[:,ny,nx] = -np.dot(mata,dpsibdx[:,ny,nx])
      
if flag_gradqbar == 1:
  qbar = np.zeros((nl,N,N))
  for nx in range(0,N):
    for ny in range(0,N):
      mata = def_radius.construct_mat(dh,gp[:,ny,nx],f0[ny,nx])
      #qbar[:,ny,nx] = -np.dot(mata,psib[:,ny,nx])/f0[ny,nx]
      qbar[:,ny,nx] = -np.dot(mata,psib[:,ny,nx])
  
  
  aux,dqbdy,dqbdx = np.gradient(qbar, Deltad)
  
# add beta effect
dqbdy = dqbdy + beta

def comp_matrices(k, l, mata, dqbdy, dqbdx, dpsibdx, dpsibdy, nu, nu4, bf):
  k2 = k**2 + l**2
  
  p2q = -mata-k2*sparse.eye(nl)
  
  diag1 = k*dqbdy - l*dqbdx
  diag1 = diag1 - (k2**2*nu + k2**3*nu4)
  diag1[-1] = diag1[-1] - k2*bf

  diag2 = l*dpsibdx - k*dpsibdy
  
  diag1 = sparse.diags(diag1,0)
  diag2 = sparse.diags(diag2,0)
  
  mat1 = diag1 + diag2.dot(p2q)
  
#  return mat1, sparse.csc_matrix(p2q)
  return mat1.toarray(), p2q.toarray()


if len(save_dat) == 0:
  save_dat = [None]*N*N
  
# # specific location
# figure middle domain: 320,320
nx1 = 20
ny1 = 20
nx2 = nx1+1
ny2 = ny1+1
nsk = 1

# nx1 = 0
# ny1 = 0
# nx2 = N
# ny2 = N
# nsk = 1

k_grid = np.zeros((N,N))
l_grid = np.zeros((N,N))
o_grid = np.zeros((N,N))

v0 = np.zeros(nl,dtype=complex)
v0 = v0 + np.random.random(nl)


ntot = len(range(nx1,nx2,nsk))*len(range(ny1,ny2,nsk))
if nx2-nx1 > 1 or ny2-ny1 > 1:
  print ('computing', ntot, 'points')
nco = 0
nco1 = 0
for nx in range(nx1,nx2,nsk):
  for ny in range(ny1,ny2,nsk):
    nco1 = nco1 + 1

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

      # periodic save
      if nco%20 == 0 and flag_savedat:
        np.save(file2, save_dat)

      nco = nco + 1
      if nx2-nx1 > 1 or ny2-ny1 > 1:
        print (nco,'; ',nco1, '/', ntot, '  ', ny + N*nx)
      # compute deformation radius
      Rd = def_radius.cal_rad(dh,gp[:,ny,nx],f0[ny,nx])
      Rd1 = Rd[1]
      
      mata = def_radius.construct_mat(dh,gp[:,ny,nx],f0[ny,nx],sparse=True)
#      l2m,m2l = def_radius.cal_transfo(dh,gp[:,ny,nx],f0[ny,nx])
           
      Rd2 = Rd1**2
      # si_k = 75
      # si_l = 37
      
      # size of the k-l window
      if doit == 1:
      # figure middle domain: 3.5
        window_s = 10.5
#        window_s = 100
      else:
        window_s = save_dat[ny + N*nx]['window_s']


      # figure middle domain: 150,75
      si_k = 50
      if (nu == 0 and nu4 == 0 and Ek == 0):
        si_l = int(si_k/2)
        l0min = 1e-3/(Rd1)
      else:
        si_l = si_k
        l0min = -window_s/(Rd1)
      
        
      k0min = -window_s/(Rd1)
      k0max = window_s/(Rd1)
      # get l<0 by symmetry
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
      
      omegai_c = np.zeros((si_l,si_k))
      eivec_c  = np.zeros((si_l,si_k,nl),dtype='complex')
      
      if flag_klmap:
        # ke = 0.1/(Rd[1]*1e3)
        # #le = 1/(Rd[1]*1e3)
        # le = 0.0 #0.1/(Rd[1]*1e3)
        for il in range(0,si_l):
          for ik in range(0,si_k):
          
            mat1_sp, mat2_sp = comp_matrices(kt[ik], lt[il], mata, dqbdy[:,ny,nx], dqbdx[:,ny,nx], dpsibdx[:,ny,nx], dpsibdy[:,ny,nx], nu, nu4, bf)
            # eival,eivec = sparse.linalg.eigs(A=sparse.linalg.inv(mat2_sp)*mat1_sp,v0=v0,k=1,which='LI')
            # omegai_c[il,ik] = -np.imag(eival[0])
            # eivec_c[il,ik,:] = np.real(eivec[:,0])
            # v0 = eivec

            eival,eivec = la.eig(mat1_sp,mat2_sp)
            idx = np.argsort(np.imag(eival))[::-1]
            omegai_c[il,ik] = np.imag(eival[idx])[0]
            eivec_c[il,ik,:] = eivec[:,idx][:,0]


              
        (ilmax,ikmax) = np.unravel_index(omegai_c[:,:].argmax(), omegai_c[:,:].shape)
        #ikmax = np.argmax(omegai_c[:,0])

        
        kmax = kt[ikmax]
        lmax = lt[ilmax]
        omax = np.real(omegai_c[ilmax,ikmax])

        eivecmax = eivec_c[ilmax,ikmax,:]
        
        fw = np.int(si_l/4)
        
        om_bdy0 = np.sum(omegai_c[-1,:]) + np.sum(omegai_c[:,0]) + np.sum(omegai_c[:,-1])
        om_bdy1 = np.sum(omegai_c[-1-fw,:]) + np.sum(omegai_c[:,+fw]) + np.sum(omegai_c[:,+fw])

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
#            print ('projection on vert modes:', np.dot(l2m[:,:],eivecmax.real))
            
        else:
          if np.sum(omegai_c[:,:]) == 0: # no instability
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
                                  'omega':omegai_c[:,:],
                                  'window_s': window_s}

      else: #if flag_klmap
        def find_omegamax(xk):
          mat1_sp, mat2_sp = comp_matrices(xk[0], xk[1], mata, dqbdy[:,ny,nx], dqbdx[:,ny,nx], dpsibdx[:,ny,nx], dpsibdy[:,ny,nx], nu, nu4, bf)

          eival,eivec = la.eig(mat1_sp,mat2_sp)
          idx = np.argsort(np.imag(eival))[::-1]

#          eival,eivec = sparse.linalg.eigs(A=sparse.linalg.inv(mat2_sp)*mat1_sp,v0=v0,k=1,which='LI')
          #TODO modify global varable v0
          #          v0 = eivec
          return -np.imag(eival[idx])[0]
        
        theta = np.arctan2(dpsibdx[0,ny,nx],-dpsibdy[0,ny,nx])
        ke = np.cos(theta)/Rd[1]
        le = np.sin(theta)/Rd[1]
        
        for ntry in [0.1,0.2,0.5,1,2,5]:
          if find_omegamax([ntry*ke,ntry*le]) !=0:
            break

        [kmax,lmax] = scipy.optimize.fmin_powell(lambda x: find_omegamax(x),[ntry*ke,ntry*le],ftol=1e-1,disp=False)
        omax = -find_omegamax([kmax,lmax])
        
        if (omax != 0.):
          k_grid[ny,nx] = kmax
          l_grid[ny,nx] = lmax
          o_grid[ny,nx] = omax
          save_dat[ny + N*nx] = {'conv':1, 
                                  'Rd':Rd,
                                  'kmax':kmax,
                                  'lmax':lmax,
                                  'omax': omax}          
  
      # 2d plot
      if flag_plot and flag_klmap:

        Rda = l_qg/(2*np.pi)
        plt.figure()
        plt.contourf(ktg*Rda,ltg*Rda,(omegai_c[:,:].real*l_qg/u_qg),20,cmap=plt.cm.hot_r)
        if (nu == 0 and nu4 == 0 and Ek == 0):
          plt.contourf(ktg*Rda,-ltg*Rda,(omegai_c[:,::-1].real*l_qg/u_qg),20,cmap=plt.cm.hot_r)
        cb = plt.colorbar(label=r'$\omega$ (t$^{-1}$)',format='%.1f')
        plt.xlabel('k (cycle/l)')
        plt.ylabel('l (cycle/l)')

        # figure middle of the domain
        #modes = 2*np.array([[48,8],[38,8],[43,2]],np.int32)

        modes = 2*np.array([[48,8]],np.int32)

        modes[0,0] = ikmax
        modes[0,1] = ilmax

#        omegai_cs = gaussian_filter(omegai_c, 3, mode='constant')
        modes = peak_local_max(omegai_c, min_distance=3)
        nmodes,naux = modes.shape

        for i in range(0,nmodes):
          ikmax2 = modes[i,1]
          ilmax2 = modes[i,0]
          plt.plot(kt[ikmax2]*Rda,lt[ilmax2]*Rda,'wo',markersize=10)
          plt.text(kt[ikmax2]*Rda,lt[ilmax2]*Rda,str(i+1),ha='center',va='center')

          
          kmax2 = kt[ikmax2]
          lmax2 = lt[ilmax2]
          omax2 = np.real(omegai_c[ilmax2,ikmax2])

          eivecmax2 = eivec_c[ilmax2,ikmax2,:]
          if flag_print:
            print ('mode:', i+1)
            print ('kmax x Rd1 = ',kmax2*Rd1)
            print ('lmax x Rd1 = ',lmax2*Rd1)
            print ('time scale = ',np.real(1/omax2/86400), ' days')
            #            print ('projection on vert modes:', np.dot(l2m[:,:],eivecmax2.real))
            
        ip = 100
        xp = np.linspace(0,2*np.pi,ip)
        psi_cos = np.repeat(np.cos(xp)[np.newaxis,:], nl, axis=0) 
        psi_sin = np.repeat(np.sin(xp)[np.newaxis,:], nl, axis=0) 
        x_tick = [0, 3.14,6.28]
        x_label = [r"$0$", r"$\pi$", r"$2\pi$"]

        for i in range(0,nmodes):
          plt.figure()
          psi = np.real(eivec_c[modes[i,0],modes[i,1],:][:,np.newaxis]*(psi_cos + 1j*psi_sin))
          plt.contourf(xp,-hh/H,psi,30,cmap=plt.cm.bwr,)
          plt.gca().set_xticks(x_tick)
          plt.gca().set_xticklabels(x_label)
          plt.xlabel("Phase")
          plt.ylabel("Depth")
          plt.title("Mode "+str(i+1))
        #plt.savefig('klo_map.png',bbox_inches='tight')
        #plt.savefig('klo_map_340_340.pdf',bbox_inches='tight')
        #plt.savefig('klo_map_340_340_30l_1.pdf',bbox_inches='tight')
        
      # 1d plot
      if flag_test2lay:
        plt.figure()
        plt.plot(kt*Rd1,omegai_c[0,:]*86400,'k')
        plt.plot(kt*Rd1,omegai_th*86400,'r--')
        plt.plot(kt*Rd1,omegai_thwob*86400,'b--')
        plt.xlim([0.0,kt[-1]*Rd1])
        
        plt.xlabel('k x Rd')
        plt.ylabel(r'$\omega (day^{-1})$')
  
if flag_savedat:
  np.save(file2, save_dat)
