#! /usr/bin/env python

import os
import re
import matplotlib
matplotlib.use('TkAgg') # set backend 
import matplotlib.pyplot as plt
#if not os.environ.has_key('DISPLAY'):
#matplotlib.use('Agg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import tkinter as Tk

from pylab import *
from struct import *
from math import pow, log10
from sys import argv, exit

class momrateclass:
    def __init__(self, fname, nt, nx, verbose=False):
        self.fname = fname
        self.nt = nt
        self.nx = nx
        self.verbose = verbose

    def readmomrate(self, xi, zi):
        self.xi = xi
        self.zi = zi
        self.V = zeros((3, self.nt), dtype=float) # xy, xx,xz
        cfname = self.fname
        fid = open(cfname, 'rb')
        nskip=(zi * self.nx + xi)
        if self.verbose:
            print("loading entry %d in %s ... " % (nskip, cfname))
        bskip= 4 * nskip * (3 + self.nt * 6) + 12
        fid.seek(bskip, 0)

        fmt="<%df" % (self.nt * 6) # > for little endian
        esiz=calcsize(fmt)

        sr = np.frombuffer(fid.read(esiz), dtype='f')
        self.V[0,:] = sr[0::6]  # xx
        self.V[1,:] = sr[5::6]  # xy
        self.V[2,:] = sr[3::6]  # xz
        fid.close()

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
         cfname="%s" % self.fname
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

if len(argv) != 2:
   print ("usage: %s directory" % argv[0])
   exit()

dirname=argv[1]


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

params = {'6.35': (354, 212), '6.45': (144, 230), '7.35': (980, 240), '8.15': (5116, 220), '8.45': (5116, 220)}
M = re.findall("\d+\.\d+", os.getcwd().split('/')[-1])[0]
(nx, nz) = params[M]
ny=4200
xskp=10
zskp = 2
nt=40000
dt=0.005
dh=0.100
#smax=0.2

srf = momrateclass("%s" % dirname, nt, nx, verbose=False)

root=Tk.Tk()
#root.title("Spontaneous Rupture Simulation")
root.title("%s" % dirname.upper() )
root.protocol("WM_DELETE_WINDOW", root.quit)

tlabel=Tk.Label(root)
tlabel.pack()
#ts=Tk.IntVar(root)
#ts.set(100)

#xs = Tk.Scale(root, from_=xskp,to=(nx-xskp), resolution=xskp,orient=Tk.HORIZONTAL)
#xs.pack(side="bottom", expand=Tk.TRUE, fill=Tk.X)
#zs=Tk.Scale(root, from_=xskp,to=(nz-xskp), resolution=xskp,orient=Tk.VERTICAL)
#zs.pack(side="right", expand=Tk.TRUE, fill=Tk.X)
xs = Tk.IntVar()
zs = Tk.IntVar()
xs.set(360)
zs.set(10)
fig = plt.figure(figsize=(6,5))

canvas=FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=Tk.TRUE)
toolbar=NavigationToolbar2Tk(canvas,root)
toolbar.pack(side="top")

def replot(event):
    fig.clf()
    xi = xs.get()
    zi = zs.get()
    srf.readmomrate(xi, zi)

    sr_xx=srf.V[0,::10]
    sr_xy=srf.V[1,::10]
    sr_xz=srf.V[2,::10]
    print ("max momrate XY: %f" % (abs(sr_xy).max()))
    #slc=srf.V[::1,yi,::1]
#fig = plt.figure(figsize=(8,8))
#   ax = fig.subplots(3, 1, sharex=True)
    ax_0 = fig.add_subplot(311)
    title("Moment rate (x=%.1f, z=%.1f)" % (xi * dh, zi * dh))
    ax_0.plot(arange(0, nt, 10) * dt, sr_xx)
    ax_0.set_ylabel("XX")
    ax_1 = fig.add_subplot(312)
    ax_1.plot(arange(0, nt, 10) * dt, sr_xy)
    ax_1.set_ylabel("XY")
    ax_2 = fig.add_subplot(313)
    ax_2.plot(arange(0, nt, 10) * dt, sr_xz)
    ax_2.set_ylabel("XZ")
    ax_2.set_xlabel("Time (s)")

#axis("image")

    canvas.draw()

slider_1=Tk.Scale(root, from_=xskp,to=(nx-xskp), resolution=xskp,
        orient='horizontal', variable=xs)
slider_1.pack(side="bottom", expand=1, fill='x')
slider_2=Tk.Scale(root, from_=zskp,to=(nz-zskp), resolution=zskp,
        orient='vertical', variable=zs)
slider_2.pack(side="left", expand=1, fill='y')
#slider.config(command=replot)
slider_1.bind('<ButtonRelease>', replot)
slider_2.bind('<ButtonRelease>', replot)

replot("<ButtonRelease>")

root.mainloop()
