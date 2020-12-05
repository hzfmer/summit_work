#! /usr/bin/env python

import matplotlib
#from pylab import *
matplotlib.use('Agg')
from matplotlib.pylab import *
import numpy as np
from os.path import isfile
from karman import karman2D_stress, karman2D_slip

dh=100.
nx=3000
nz=180

xis=2600
#xis=300
zis=100

# shallower source: scenario Z
# zis = 50

zi1=170
xi0=200
xi1=2700

#taper on mu_d for
tpw_h=100 #horizontals
tpw_b=50  #bottom
tpw_s=30  #surface

# A & B & G 
mu_stat=0.55
mu_dyn=0.325
# C & D & H
#mu_stat=0.65
#mu_dyn=0.3
# E & F & I
#mu_stat = 0.55
#mu_dyn = 0.315

# J & K & L
#mu_stat = 0.55
#mu_dyn = 0.358

# M & N & O
#mu_stat = 0.55
#mu_dyn = 0.348

#mu_stat = 0.6
#mu_dyn = 0.3

# E & F & I (discarded)
#ratio = 0.41 
# A & B & C & D & G & H & Z
ratio = 0.4
# J (discarded)
# ratio = 0.44
axm=47.0e3
azm=9.0e3
H=0.50

def compute_vertical_stress(rho, dh, G):
   nz,nx=rho.shape
   zz=np.zeros((nz,nx), dtype=float)
 
   zz[0,:] = -dh/2. * rho[0,:] * G

   for k in range(1, nz):
      zz[k,:] = zz[k-1,:] - dh * rho[k,:] * G
   
   return zz

fl=np.loadtxt("fault_3000.dat")
rho=fl[:,5].reshape(nz,nx)
vstmp=fl[:,4].reshape(nz,nx)
vs=np.zeros((nz,nx), dtype=float)
for k in range(nz):
   depth = k * dh
   if depth < 4000.:
      zfact = 0.7
   elif depth >= 4000. and depth <= 6000.:
      zfact = 0.7 + (depth - 4000.)/2000.* 0.3;
   else:
      zfact = 1.
   for l in range(nx):
      vs[k,l] = max(500., vstmp[k,l] * zfact)
rho_w=np.ones((nz,nx), dtype=float) * 1000.
del fl
dpt=np.arange(0, nz*dh, dh)

zz=compute_vertical_stress(rho, dh, 9.81)
pf=-compute_vertical_stress(rho_w, dh, 9.81)

yy=zz+pf
xy=(zz+pf) * (-ratio)

mu_s=np.ones((nz,nx), dtype=float) * mu_stat
mu_d=np.ones((nz,nx), dtype=float) * mu_dyn
dmax=np.ones((nz,nx), dtype=float) * 10.

#mu_s[0:zi1,xi0:xi1]=mu_stat
#mu_d[0:zi1,xi0:xi1]=mu_dyn
dmax[0:zi1,xi0:xi1]=0.3

coh=np.zeros((nz,), dtype=float)
tnuc=np.zeros((nz,nx), dtype=float)

if not isfile("karman.npy"):
   Fc=karman2D_slip(nx, nz, dh, axm, azm, H)
   if abs(np.real(Fc).min()) < np.real(Fc).max():
      Fc=-Fc
   np.save("karman", Fc)
else:
   Fc=np.load("karman.npy")
 
F=np.real(Fc)
F/=abs(F).max()
print F.min(), F.max()
#imshow(F)

#apply random field to mu_s, mu_d:
for k in range(zi1):
   for l in range(xi0,xi1):
       mu_s[k,l] -= (mu_s[k,l] - ratio) * F[k,l]
       mu_d[k,l] -= (mu_s[k,l] - ratio) * F[k,l]

#taper off mu_d along edges of fault:
#bottom
for l in range(xi0, xi1):
   for k in range(zi1-tpw_b, zi1):
      #mu_d[k,l]=mu_d[k,l] + float(k-zi1+tpw_b) / float(tpw_b) * (mu_s[k,l]-mu_d[k,l])
      dmax[k,l]=0.3 + float(k-zi1+tpw_b) / float(tpw_b) * 9.7

#surface
#for l in range(xi0, xi1):
#   for k in range(0, tpw_s):
#      mu_d[k,l]=mu_d[k,l]+ float(tpw_s-k) / float(tpw_s) * (mu_s[k,l]-mu_d[k,l])

#left
for l in range(xi0, xi0+tpw_h):
   for k in range(0, zi1):
      #mu_d[k,l]=max(mu_d[k,l], mu_s[k,l] - float(l-xi0) / float(tpw_h) * (mu_s[k,l]-mu_d[k,l]))
      dmax[k,l]=max(dmax[k,l], 10. - float(l-xi0) / float(tpw_h) * 9.7)

#right
for l in range(xi1-tpw_h, xi1):
   for k in range(0, zi1):
      #mu_d[k,l]=max(mu_d[k,l], mu_d[k,l] + float(l-xi1+tpw_h) / float(tpw_h) * (mu_s[k,l]-mu_d[k,l]))
      dmax[k,l]=max(dmax[k,l], 0.3 + float(l-xi1+tpw_h) / float(tpw_h) * 9.7)

     

sd=np.zeros((nz,nx), dtype=float)
se=np.zeros((nz,nx), dtype=float)

for k in range(nz):
   #if dpt[k] <= 2400.:
   #   coh[k]=0.000425e6 * (2400. - dpt[k])
   if k < tpw_s:
      coh[k]=float(tpw_s - k) / float(tpw_s) * 1.e6
   r_crit=7500.
   for l in range(nx):
      rad=sqrt(pow(xis-l,2.)+pow(zis-k,2.)) * dh

      if rad < r_crit:
          tnuc[k,l]= rad / (0.7 *vs[k,l]) + 0.081*r_crit / (0.7*vs[k,l]) * \
             ( 1. / (1. - pow(rad/r_crit, 2.)) - 1.)
      else:
          tnuc[k,l]=1.e9

      sd[k,l]=abs(xy[k,l]) - abs(yy[k,l]*mu_d[k,l])
      se[k,l]=abs(yy[k,l]*mu_s[k,l]) - abs(xy[k,l])

Sf=se/sd
print Sf[0:zi1,xi0:xi1].mean()
print "average stress drop: %e" % sd[0:zi1-tpw_b,xi0:xi1].mean()
imshow(tnuc)
#show()
savefig('tnuc.pdf')
#if 0:
#  figure()
#  plot(coh, dpt)
#  ylim(dpt.max(), dpt.min())
#  #show()
#
#if 1:
#  figure()
#  imshow(Sf, vmin=0.4, vmax=2.0)
#  colorbar()
#
#if 1:
#  figure()
#  plot (-yy[:,500], dpt, label='yy')
#  plot (xy[:,500], dpt, label='xy')
#  plot (-yy[:,500]*mu_s[:,500]+coh, dpt, label='static shear strength')
#  plot (-yy[:,500]*mu_d[:,500]+coh, dpt, label='dynamic shear strength')
#  axis((0, (-yy).max(), dpt.max(), dpt.min()))
#  legend(loc=1)
#  #show()
#
#if 0:
#  figure()
#  #bins=arange(0, 5.e7, 5.e5)
#  #hist(sd[80,xi0+tpw_h:xi1-tpw_h], bins)
#  bins=arange(0.4, .6, 1.e-3)
#  hist(mu_s[:,xi0+tpw_h:xi1-tpw_h].flatten(), bins)
#  show()
#
#if 1:
#  figure(figsize=(15,8))
#
#  x=linspace(-nx/2 * dh, nx/2 * dh, nx) / 1.e3
#  z=dpt / 1.e3
#
#  subplot(231)
#  #imshow(mu_s, vmin=0, vmax=1.)
#  contourf(x, z, mu_s, arange(0.3, 0.7, 0.01), extend="both")
#  axis("image")
#  axis((x.min(), x.max(), z.max(), 0))
#  #clim(0, 1)
#  title('$\\mu_s$')
#  colorbar(orientation="horizontal")
#
#  subplot(232)
#  #imshow(mu_d, vmin=0, vmax=1)
#  contourf(x, z, mu_d, arange(0.2, 0.6, 0.01), extend="both")
#  axis("image")
#  axis((x.min(), x.max(), z.max(), 0))
#  #clim(0, 1)
#  title('$\\mu_d$')
#  colorbar(orientation="horizontal")
#
#  subplot(233)
#  #imshow(mu_d, vmin=0, vmax=1)
#  contourf(x, z, sd, arange(0., 1.e7, 1.e5), extend="both")
#  axis("image")
#  axis((x.min(), x.max(), z.max(), 0))
#  #clim(0, 1)
#  title('Stress Drop')
#  colorbar(orientation="horizontal")
#
#
#  subplot(234)
#  clim=arange(0, 1., 0.01)
#  contourf(x, z, tnuc, clim)
#  axis("image")
#  axis((x.min(), x.max(), z.max(), 0))
#  title('$T_\mathrm{nuc}$')
#  colorbar(orientation="horizontal", ticks=(0,1))
#
#  subplot(235)
#  clim=arange(0, 10.1, 0.01)
#  contourf(x, z, dmax, clim)
#  axis("image")
#  axis((x.min(), x.max(), z.max(), 0))
#  title('$D_\mathrm{max}$')
#  colorbar(orientation="horizontal", ticks=(0,1))
#
#
#  show()

fid=open("plast_M7.8_z0_cvm_D_lvz2_A.dat", "w")
yi=1172
x0=0
yz=0.
for k in range(nz):
   for l in range(20, nx-20):
     fid.write("%5d %5d %5d %e %e %e %e %e %e %e %e\n" % \
       (l+1+x0, yi+1, k+1, -yy[k,l], xy[k,l], dmax[k,l], yz, mu_s[k,l], mu_d[k,l], coh[k], tnuc[k,l]))

fid.close()
