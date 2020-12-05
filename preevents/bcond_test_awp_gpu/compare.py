#! /usr/bin/env python

import pylab as pl

vel1=pl.loadtxt("output_sfc/SX_0_0032_0032_0001.dat")
dt=0.0025

#vel1=-pl.cumsum(vel1)*d1
nt=len(vel1)
t=pl.linspace(dt, nt*dt, nt)

vel2=pl.loadtxt("../bcond_test/output_sfc/SX0032_0032_0064.dat")

pl.plot(t, vel1, label="AWP GPU", linewidth=3, color='gray')
pl.plot(t, vel2, label="AWP CPU", color='blue')
pl.legend()
pl.show()
