#! /usr/bin/env python

import matplotlib
matplotlib.use('TkAgg') # set backend 
#matplotlib.use('Agg') # set backend 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg

from pylab import *
from struct import *
from math import pow, log10
from sys import argv, exit

class srfvelclass:
   def __init__(self, fname, nx, ny, nz, dh, tskip, wstep, verbose=False):
      self.fname=fname
      self.nx=nx
      self.ny=ny
      self.nz=nz
      self.dh=dh
      self.tskip=tskip
      self.wstep=wstep
      self.xoff=0.
      self.yoff=0.
      self.verbose=verbose
   
   def readsurf(self, ti):
      self.ti=ti
      askip=self.wstep * self.tskip;

      self.V=zeros( (self.nz,self.ny,self.nx), dtype=float)

      if self.tskip==0:
         cfname=self.fname
         fid=open(cfname, "r")

      elif self.wstep==0:
         cfname=self.fname
         fid=open(cfname, "r")
         nskip=(ti / self.tskip) - 1
         if self.verbose:
            print("loading entry %d in %s ... " % (nskip, cfname)) 
         bskip=calcsize("<%df" % nskip) * self.nx * self.ny * self.nz
         fid.seek(bskip, 0)
      else:
         fnum=ceil(ti/float(askip)) * askip
         cfname="%s%07d" % (self.fname, fnum )
         nskip=float(ti-fnum+askip-1 )/float(self.tskip)
         if self.verbose:
            print ("loading entry %d in %s ... " % (nskip, cfname)) 
         bskip=calcsize("<%df" % nskip) * self.nx * self.ny * self.nz
         fid=open(cfname, "r")
         fid.seek(bskip, 0)

      fmt="<%df" % self.nx # > for little endian
      esiz=calcsize(fmt)

      for m in range(self.nz):
        for n in range(self.ny):
          sr=fid.read(esiz)
          self.V[m,n,:]=unpack(fmt, sr)
   
      fid.close()

if len(argv) != 3:
   print ("usage: %s directory component(X/Y/Z)" % argv[0])
   exit()

filename=argv[1] + "/S" + argv[2] + "_0_"

import Tkinter as Tk

def readslip(nx, nz, fname):
   npt=nx*nz

   fmt="<%df" % npt
   siz=calcsize(fmt)

   fid=open(fname, "r")
   asc=fid.read(siz)
   fid.close()

   dat=unpack(fmt, asc)
   slip=reshape(dat, (nz,nx))
   return slip

def readrate19(nx, ny, nz, fbase, timestep, tskp, wstep, var):
   params=["rateu_u", "rateu_w", "slipu_u", "slipu_w", "tru1", "trw1", "trv1", \
           "uplus", "vplus", "wplusi", "umin", "vmin", "wmin", "duplus", "dvplus", \
           "dwplus", "dumin", "dvmin", "dwmin"]
   idx=params.index(var)

   askip=tskp * wstep
   fnum=ceil(timestep/float(askip)) * askip
   npt=nx*ny*nz

   fmt="<%df" % npt
   siz=calcsize(fmt)

   fname="%s%07d" % (fbase, fnum)

   nskip=float(timestep-fnum+askip)/float(tskp) - 1
   print ("loading entry %f in %s ... " % (nskip, fname))

   fid=open(fname, "r")
   bskip=calcsize("<1f") * (nx*ny*nz*nskip*19 + nx*ny*nz*idx)
   fid.seek(bskip, 0)

   asc=fid.read(siz)
   fid.close()

   dat=unpack(fmt, asc)
   rate=reshape(dat, (nz,ny,nx))
   return rate

def readrate2(nx, ny, nz, fname, timestep, tskp):
   npt=nx*ny*nz

   fmt="%df" % npt
   siz=calcsize(fmt)

   if timestep % tskp != 0:
      print ("warning: timestep not divisible by tskp")
   nskip=timestep / tskp - 1
   print ("loading entry %f in %s ... " % (nskip, fname))

   fid=open(fname, "r")
   bskip=calcsize("<1f") * nx * ny * nz * nskip
   fid.seek(bskip, 0)

   asc=fid.read(siz)
   fid.close()

   dat=unpack(fmt, asc)
   rate=reshape(dat, (nz,ny,nx))
   rate=flipud(rate)
   return rate

#nproc=int(dirname.split('.')[1].split('/')[0])

nx=64 #/ nproc
ny=64
nz=128
tskp=10
nt=4000
wstep=20
dt=0.0025
dh=1.00
smax=2.
zi=54
xi=32
yi=32

srf=srfvelclass("%s" % filename, nx, ny, nz, dh, tskp, wstep, verbose=False)

x=linspace(0., dh*(nx-1), nx)
y=linspace(0., dh*(ny-1), ny)
z=linspace(0., dh*(nz-1), nz)
#z=flipud(z) #used only for SeismIO
x=x-x.mean()
y=y-y.mean()

root=Tk.Tk()
#root.title("Spontaneous Rupture Simulation")
root.title("%s" % filename.upper() )
root.protocol("WM_DELETE_WINDOW", root.quit)

tlabel=Tk.Label(root)
tlabel.pack()
ts=Tk.IntVar(root)
ts.set(1000)

fig=figure(figsize=(8,8))

canvas=FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=Tk.TRUE)
toolbar=NavigationToolbar2TkAgg(canvas,root)
toolbar.pack(side="top")

def replot(event):
   fig.clf()
   t=ts.get()
   fig.add_subplot(111)
   srf.readsurf(t)

   print ("min, max X: %f, %f" % (srf.V.min(), srf.V.max()))
   #slc=srf.V[:,yi,:]
   slc=srf.V[:,:,xi]
   #slc=srf.V[zi,:,:]
   smax=abs(slc).max()
   #imshow(slc, vmin=-smax, vmax=smax)
   cnt=linspace(-smax, smax, 50)
   contourf(y, z, slc, cnt, extend="both")
   xlabel("Y")
   ylabel("Z")
   cab = colorbar(orientation="horizontal", shrink=0.5, ticks=(-smax, 0., smax))
   cab.set_label("$m/s$", fontsize=20)
   title("V%s (%5.3f sec)" % (argv[2], t*dt))
   axis((y.min(), y.max(), z.max(), 0.))
   axis("image")

   savefig("snap_wwl_12_7040.pdf")
   canvas.draw()

slider=Tk.Scale(root, from_=tskp,to=(nt-tskp), resolution=tskp,orient=Tk.HORIZONTAL, variable=ts)
slider.pack(side="left", expand=Tk.TRUE, fill=Tk.X)
#slider.config(command=replot)
slider.bind('<ButtonRelease>', replot)

replot("<ButtonRelease>")

root.mainloop()
