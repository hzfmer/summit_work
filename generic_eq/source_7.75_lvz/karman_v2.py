#! /usr/bin/env python

#2D Karman ACF based on a Matlab script by J. Ripperger
#nx, nz: grid dimension (in samples)
#samp: sampling distance
#axm, azm: along-strike and along-dip autocorrelation length (same unit as samp)
#H: Hurst exponent

from pylab import *

def karman(nx, nz, samp, axm, azm, H):
   knyx=1./(2*samp)
   knyz=1./(2*samp)

   if nx % 2 == 0:
      rmx=nx/2
   else:
      rmx=(nx-1)/2
   if nz % 2 == 0:
      rmz=nz/2
   else:
      rmz=(nz-1)/2

   kx=linspace(-knyx, knyx, 2*rmx+1)
   kz=linspace(-knyz, knyz, 2*rmz+1)

   PS=zeros((rmz+1,2*rmx+1), dtype=float)

   for k in range(rmz+1):
     for l in range(2*rmx+1):
        kr = sqrt(pow(kz[k]*azm,2.) + pow(kx[l]*axm,2.)) 
        PS[k,l] = axm * azm / (1. + pow(pow(kr, 2.), H+1.))

   PS /= PS.max()

   AM = sqrt(PS)
   PH = 2*pi*random(PS.shape)
   FF = AM*cos(PH) + 1j * AM * sin(PH)
   FF2=append(FF.copy(), conj(fliplr(flipud(FF.copy()[0:rmz,:]))), axis=0)

   for k in range(rmx):
      FF2[rmz,-k+2*rmx] = conj(FF2[rmz,k])

   F=ifftn(ifftshift(FF2)).real

   return(F)

"""
nx=500
nz=500

samp = 1.
axm=500.
azm=100.
H = 0.75
F=karman(nx, nz, samp, axm, azm, H)

contourf(real(F))
axis("image")
colorbar()
show()
"""
