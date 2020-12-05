#! /usr/bin/env python

import os
import matplotlib
matplotlib.use('TkAgg') # set backend 
#if not os.environ.has_key('DISPLAY'):
#matplotlib.use('Agg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg

from pylab import *
from struct import *
from math import pow, log10
from sys import argv, exit
import numpy as np

if len(argv) != 2:
   print "usage: %s directory" % argv[0]
   exit()

fname=argv[1]

import Tkinter as Tk

#nproc=int(dirname.split('.')[1].split('/')[0])

nx=1400 
ny=1400
nz=600
dh=0.020
nvar = 3
#smax=0.2


root=Tk.Tk()
#root.title("Spontaneous Rupture Simulation")
root.title(" Slice along X direction")
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

Vtmp = np.fromfile(fname,'f').reshape(nz, ny, nx, nvar)
def replot(event):
   fig.clf()
   t=int(ts.get())
   fig.add_subplot(111)
   V = Vtmp[:,:,t,2]
   print "min, max X: %f, %f" % (V.min(), V.max())
   slc=V[::1,:]
   smax=abs(slc).max()
   imshow(slc, vmin=0, vmax=smax)

   colorbar(orientation="horizontal")
   title("X = %d" % t)
   axis("image")

#   savefig("snap_SX5000.pdf")
   canvas.draw()

slider=Tk.Scale(root, from_= 0,to=nx, resolution=1,orient=Tk.HORIZONTAL, variable=ts)
slider.pack(side="left", expand=Tk.TRUE, fill=Tk.X)
#slider.config(command=replot)
slider.bind('<ButtonRelease>', replot)

replot("<ButtonRelease>")

root.mainloop()
