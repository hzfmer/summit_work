#! /usr/bin/env python

from pylab import *

def karman2D_slip(nx, nz, samp, axm, azm, H):
   PS=zeros((nz,nx), dtype=float)
   PH=zeros((nz,nx), dtype=float)
   FF=zeros((nz,nx), dtype=complex)
   kr=zeros((nz,nx), dtype=float)

   knyx=1./(2*samp)
   knyz=1./(2*samp)
   rmx=nx/2
   rmz=nz/2

   kx=linspace(-knyx, knyx, nx)
   kz=linspace(-knyz, knyz, nz)

   for k in range(nz):
     for l in range(nx):
        kr[k,l] = sqrt(pow(kz[k]*azm,2.) + pow(kx[l]*axm,2.)) 
        PS[k,l] = axm * azm / (1. + pow(pow(kr[k,l], 2.), H+1.))
        PH[k,l] = 2*pi*rand()

   for k in range(rmz):
     for l in range(nx):
        am = sqrt(PS[k,l])
        FF[k,l] = am*cos(PH[k,l]) + 1j * am * sin(PH[k,l])

   for k in range(rmz+1,nz):
     for l in range(nx):
        FF[k,l] = conj(FF[nz-k-1,nx-l-1])


   F=ifftn(ifftshift(FF))

   return(F)

def karman2D_stress(nx, nz, samp, axm, azm, H):
   PS=zeros((nz,nx), dtype=float)
   PH=zeros((nz,nx), dtype=float)
   FF=zeros((nz,nx), dtype=complex)
   kr=zeros((nz,nx), dtype=float)

   knyx=1./(2*samp)
   knyz=1./(2*samp)
   rmx=nx/2
   rmz=nz/2

   kx=linspace(-knyx, knyx, nx)
   kz=linspace(-knyz, knyz, nz)

   for k in range(nz):
     for l in range(nx):
        kr[k,l] = sqrt(pow(kz[k]*azm,2.) + pow(kx[l]*axm,2.)) 
        PS[k,l] = axm * azm *(pow(kx[l], 2.) + pow(kz[k], 2.)) \
           / (1. + pow(pow(kr[k,l], 2.), H+1.))
        PH[k,l] = 2*pi*rand()

   for k in range(rmz):
     for l in range(nx):
        am = sqrt(PS[k,l])
        FF[k,l] = am*cos(PH[k,l]) + 1j * am * sin(PH[k,l])

   for k in range(rmz+1,nz):
     for l in range(nx):
        FF[k,l] = conj(FF[nz-k-1,nx-l-1])


   F=ifftn(ifftshift(FF))

   return(F)



"""
nx=201
nz=101

samp = 100.
axm=10000.
azm=500.
H = 0.75
F=karman2D(nx, nz, samp, axm, azm, H)

contourf(real(F))
axis("image")
colorbar()
show()
"""
