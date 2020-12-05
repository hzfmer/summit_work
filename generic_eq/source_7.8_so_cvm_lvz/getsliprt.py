#! /usr/bin/env python

from pylab import *
from struct import *
from sys import stdout

nt=5000
nflt=1560
fname="momrate_LaHabra.zf.100m.bin"
dt=0.001

def getmom2(dat):
   mij_squared=pow(dat[0], 2.) + pow(dat[1], 2.) + pow(dat[2], 2.) + \
      2.*pow(dat[3], 2.) + 2.*pow(dat[4], 2.) + 2.*pow(dat[5], 2.)
   mom=1./sqrt(2.) * sqrt(mij_squared)
   return mom

def getmulmom2(dat,nt):
   mom=zeros((nt,), dtype=float)
   for k in range(nt):
      mij_squared=pow(dat[k,0], 2.) + pow(dat[k,1], 2.) + pow(dat[k,2], 2.) + \
           2.*pow(dat[k,3], 2.) + 2.*pow(dat[k,4], 2.) + 2.*pow(dat[k,5], 2.)
      mom[k]=1./sqrt(2.) * sqrt(mij_squared)
   return mom

def subfault_mom_rt(dat, dt):
   nt,nel=dat.shape   
   #for n in range(nt):
   #   mom=getmom(dat[n,:])
   #   fltmom+=mom
   #   if mom > 0 and rt == 0:
   #      rt=n
   #smom=sum(dat, axis=0) * dt
   #fltmom=getmom(smom)
   #mom=getmulmom2(dat,nt)
   #fltmom=sum(mom*dt)
   idx=where(abs(dat).max(axis=0) > 0.)
   if len(idx[0]) > 0:
      rt=idx[0][0]
   else:
      rt=-999
   return fltmom,rt

def ff_moment_rt(fname, nflt, nt, dt):
   fid=open(fname, "r")
   pos=zeros((nflt,3), dtype=int)
   mom=zeros((nflt,), dtype=float)
   rt=zeros((nflt,), dtype=int)

   hfmt="3i"
   dfmt="%df" % (6*nt)

   for k in range(nflt):
      stdout.write("\rextracting subfault %d out of %d" % (k+1, nflt))
      stdout.flush()
      hasc=fid.read(calcsize(hfmt))
      dasc=fid.read(calcsize(dfmt)) 

      pos[k,:]=unpack(hfmt, hasc)
      buf=array(unpack(dfmt, dasc)).reshape(nt, 6)
      mom[k],rt[k]=subfault_mom_rt(buf, dt)

   fid.close()
   return pos, mom, rt

def write_moment_rt(pos,mom,rt,nflt,dt):
   fid=open("momrate_la_habra_zf_100m.asc", "w")
   for k in range(nflt):
      rta=rt[k]*dt
      fid.write("%d %d %d %e %f\n" % (pos[k,0], pos[k,1], pos[k,2], mom[k], rta))
   fid.close()

pos,mom,rt=ff_moment_rt(fname, nflt, nt, dt)
write_moment_rt(pos,mom,rt,nflt,dt)

def mom2mw(m0):
   mw=2./3. * (log10(m0) -9.1)
   return mw

m0=sum(mom)

print "Total Moment: %e Nm (Mw = %4.1f)" % (m0, mom2mw(m0))


#imshow(mom.reshape(150,260))
#axis((0, 260, 0, 150))
#colorbar()
#show()

