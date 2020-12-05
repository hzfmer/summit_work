#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import os
import matplotlib
matplotlib.use('TkAgg') # set backend 
#if not os.environ.has_key('DISPLAY'):
#matplotlib.use('Agg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

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
         cfname="%s_0_%07d" % (self.fname, fnum )
         nskip=float(ti-fnum+askip-1 )/float(self.tskip)
         print(ti, fnum, askip, self.tskip)
         if self.verbose:
            print ("loading entry %d in %s ... " % (nskip, cfname)) 
         bskip=calcsize("<%df" % nskip) * self.nx * self.ny * self.nz
         fid=open(cfname, "rb")
         fid.seek(bskip, 0)

      fmt="<%df" % self.nx # > for little endian
      esiz=calcsize(fmt)

      for m in range(self.nz):
        for n in range(self.ny):
          sr=fid.read(esiz)
          self.V[m,n,:]=unpack(fmt, sr)
   
      fid.close()

if len(argv) != 2:
   print ("usage: %s directory" % argv[0])
   exit()

dirname=argv[1]

import tkinter as Tk

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

   fname="%s_0_%07d" % (fbase, fnum)

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

nx=3456 #/ nproc
ny=3456
nz=1
tskp=20
nt=20000
#wstep=100
wstep=200
dt=0.0005
dh=0.008
#smax=0.2
zi=0

srf=srfvelclass("%s/SX" % dirname, nx, ny, nz, dh, tskp, wstep, verbose=True)

root=Tk.Tk()
#root.title("Spontaneous Rupture Simulation")
root.title("%s" % dirname.upper() )
root.protocol("WM_DELETE_WINDOW", root.quit)

tlabel=Tk.Label(root)
tlabel.pack()
ts=Tk.IntVar(root)
ts.set(1000)

fig=figure(figsize=(8,8))

canvas=FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=Tk.TRUE)
toolbar=NavigationToolbar2Tk(canvas,root)
toolbar.pack(side="top")

def replot(event):
   fig.clf()
   t=ts.get()
   fig.add_subplot(111)
   srf.readsurf(t)

   print ("min, max X: %f, %f" % (srf.V.min(), srf.V.max()))
   slc=srf.V[zi,::-8,::8]
   #slc=srf.V[::1,yi,::1]
   smax=abs(slc).max()
   imshow(slc, vmin=-smax, vmax=smax)

   colorbar(orientation="horizontal")
   title("X (%5.3f sec)" % (t*dt))
   axis("image")

   # savefig("snap_SX5000.pdf")
   canvas.draw()

slider=Tk.Scale(root, from_=tskp,to=(nt-tskp), resolution=tskp,orient=Tk.HORIZONTAL, variable=ts)
slider.pack(side="left", expand=Tk.TRUE, fill=Tk.X)
#slider.config(command=replot)
slider.bind('<ButtonRelease>', replot)

replot("<ButtonRelease>")

root.mainloop()
