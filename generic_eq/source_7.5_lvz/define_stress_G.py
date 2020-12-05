#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np
from os.path import isfile
from karman_v2 import karman

dh=100.
nx=1920
nz=200

xis=1400
zis=100

zi1=170
xi0=400
xi1=1500

#taper on mu_d for
tpw_h=100 #horizontals
tpw_b=50  #bottom
tpw_s=30  #surface

mu_stat=0.60
mu_dyn=0.30

axm=31.26e3
azm=7.16e3
H=0.75

def compute_vertical_stress(rho, dh, G):
   nz=len(rho)
   zz=zeros((nz,), dtype=float)
 
   zz[0] = -dh/2. * rho[0] * G

   for k in range(1, nz):
      zz[k] = zz[k-1] - dh * rho[k] * G
   return zz

par=loadtxt("../velmods/mojave_100m.dat")
#dpt,vp,vs,rho,qp,qs
dpt=par[0:nz,0]
vs=par[0:nz,2]
rho=par[0:nz,3]
nz=len(dpt)

del par

rho_w=ones((nz,), dtype=float) * 1000.

zz=compute_vertical_stress(rho, dh, 9.81)
pf=-compute_vertical_stress(rho_w, dh, 9.81)

yy=zz+pf
xy=(zz+pf) * (-0.4)

mu_s=ones((nz,nx), dtype=float) * 1.e5
mu_d=ones((nz,nx), dtype=float) * 1.e5

mu_s[0:zi1,xi0:xi1]=mu_stat
mu_d[0:zi1,xi0:xi1]=mu_dyn


coh=zeros((nz,), dtype=float)
tnuc=zeros((nz,nx), dtype=float)

if not isfile("karman.npy"):
   Fc=karman(nx, nz, dh, axm, azm, H)
   if abs(real(Fc).min()) < real(Fc).max():
      Fc=-Fc
   np.save("karman", Fc)
else:
   Fc=np.load("karman.npy")
 
F=real(Fc)
F/=abs(F).max()
print F.min(), F.max()
#imshow(F)

#apply random field to mu_s, mu_d:
for k in range(zi1):
   for l in range(xi0,xi1):
       mu_s[k,l] -= (mu_s[k,l] - 0.4) * F[k,l]
       mu_d[k,l] -= (mu_s[k,l] - 0.4) * F[k,l]

#taper off mu_d along edges of fault:
#bottom
for l in range(xi0, xi1):
   for k in range(zi1-tpw_b, zi1):
      mu_d[k,l]=mu_d[k,l] + float(k-zi1+tpw_b) / float(tpw_b) * (mu_s[k,l]-mu_d[k,l])

#surface
#for l in range(xi0, xi1):
#   for k in range(0, tpw_s):
#      mu_d[k,l]=mu_d[k,l]+ float(tpw_s-k) / float(tpw_s) * (mu_s[k,l]-mu_d[k,l])

#left
for l in range(xi0, xi0+tpw_h):
   for k in range(0, zi1):
      mu_d[k,l]=max(mu_d[k,l], mu_s[k,l] - float(l-xi0) / float(tpw_h) * (mu_s[k,l]-mu_d[k,l]))

#right
for l in range(xi1-tpw_h, xi1):
   for k in range(0, zi1):
      mu_d[k,l]=max(mu_d[k,l], mu_d[k,l] + float(l-xi1+tpw_h) / float(tpw_h) * (mu_s[k,l]-mu_d[k,l]))

     

dmax=0.94
sd=zeros((nz,nx), dtype=float)
se=zeros((nz,nx), dtype=float)

for k in range(nz):
   #if dpt[k] <= 2400.:
   #   coh[k]=0.000425e6 * (2400. - dpt[k])
   if k < tpw_s:
      coh[k]=float(tpw_s - k) / float(tpw_s) * 1.e6
   r_crit=10000.
   for l in range(nx):
      rad=sqrt(pow(xis-l,2.)+pow(zis-k,2.)) * dh

      if rad < r_crit:
          tnuc[k,l]= rad / (0.7 *vs[k]) + 0.081*r_crit / (0.7*vs[k]) * \
             ( 1. / (1. - pow(rad/r_crit, 2.)) - 1.)
      else:
          tnuc[k,l]=1.e9

      sd[k,l]=abs(xy[k]) - abs(yy[k]*mu_d[k,l])
      se[k,l]=abs(yy[k]*mu_s[k,l]) - abs(xy[k])

Sf=se/sd
print Sf[0:zi1,xi0:xi1].mean()
print "average stress drop: %e" % sd[0:zi1-tpw_b,xi0:xi1].mean()

if 0:
  figure()
  plot(coh, dpt)
  ylim(dpt.max(), dpt.min())
  #show()

if 1:
  figure()
  imshow(Sf, vmin=0.4, vmax=2.0)
  colorbar()

if 1:
  figure()
  plot (-yy, dpt, label='yy')
  plot (xy, dpt, label='xy')
  plot (-yy*mu_s[:,500]+coh, dpt, label='static shear strength')
  plot (-yy*mu_d[:,500]+coh, dpt, label='dynamic shear strength')
  axis((0, (-yy).max(), dpt.max(), dpt.min()))
  legend(loc=1)
  #show()

if 0:
  figure()
  #bins=arange(0, 5.e7, 5.e5)
  #hist(sd[80,xi0+tpw_h:xi1-tpw_h], bins)
  bins=arange(0.4, .6, 1.e-3)
  hist(mu_s[:,xi0+tpw_h:xi1-tpw_h].flatten(), bins)
  show()

if 1:
  figure(figsize=(15,8))

  x=linspace(-nx/2 * dh, nx/2 * dh, nx) / 1.e3
  z=dpt / 1.e3

  subplot(221)
  #imshow(mu_s, vmin=0, vmax=1.)
  contourf(x, z, mu_s, arange(0.3, 0.7, 0.01), extend="both")
  axis("image")
  axis((x.min(), x.max(), z.max(), 0))
  #clim(0, 1)
  title('$\\mu_s$')
  colorbar(orientation="horizontal")

  subplot(222)
  #imshow(mu_d, vmin=0, vmax=1)
  contourf(x, z, mu_d, arange(0.2, 0.6, 0.01), extend="both")
  axis("image")
  axis((x.min(), x.max(), z.max(), 0))
  #clim(0, 1)
  title('$\\mu_d$')
  colorbar(orientation="horizontal")

  subplot(223)
  #imshow(mu_d, vmin=0, vmax=1)
  contourf(x, z, sd, arange(0., 1.e7, 1.e5), extend="both")
  axis("image")
  axis((x.min(), x.max(), z.max(), 0))
  #clim(0, 1)
  title('Stress Drop')
  colorbar(orientation="horizontal")


  subplot(224)
  clim=arange(0, 1., 0.01)
  contourf(x, z, tnuc, clim)
  axis("image")
  axis((x.min(), x.max(), z.max(), 0))
  title('$T_\mathrm{nuc}$')
  colorbar(orientation="horizontal", ticks=(0,1))

  show()

fid=open("plast_M7.5_z0_G_lvz2.dat", "w")
yi=250
#fid=open("plast_M7.5_large_z0_G.dat", "w")
#yi=800
yz=0.
for k in range(nz-20):
   for l in range(20,nx-20):
     fid.write("%5d %5d %5d %e %e %e %e %e %e %e %e\n" % \
       (l+1, yi+1, k+1, -yy[k], xy[k], dmax, yz, mu_s[k,l], mu_d[k,l], coh[k], tnuc[k,l]))

fid.close()
