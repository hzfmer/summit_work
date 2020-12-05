#! /usr/bin/env python

import pylab as pl
import struct

def ricker(f, length=0.512, dt=0.001):
    t = pl.linspace(-length/2, (length-dt)/2, length/dt)
    y = (1.0 - 2.0*(pl.pi**2)*(f**2)*(t**2)) * pl.exp(-(pl.pi**2)*(f**2)*(t**2))
    return t, y

def vel2momrate(vel, mu, dh, dt):
    nt=len(vel)
    mr=pl.zeros((nt,), dtype=float)
    for n in range(nt):
       mr[n] = 2.*mu*pow(dh,2.)*vel[n]
    return mr

dt=0.00025
Vs=500.
rho=2000.
dh=1.
f=10.
mu=pow(Vs,2.)*rho

Vp=1500.

nx=64
ny=64
nz=64

pmlw=10

length=0.25

t,velx=ricker(f, length=0.25, dt=dt)
nt=len(velx)

mrx=vel2momrate(velx, mu, dh, dt)

mr=pl.zeros((nt,6), dtype=float)
mr[:,3] = mrx

mrasc=struct.pack("%df" % (6*nt), *mr.flatten())

fid=open("plane_wave3.bin_0", "w")

cnt=0
for l in range(0, ny):
   for k in range(0, nx):
      headasc=struct.pack("3i", k+1, l+1, nz-pmlw)
      fid.write(headasc)
      fid.write(mrasc)
      cnt+=1
fid.close()

print("nst=%d, nsrc=%d" % (nt, cnt))
