#! /usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
import Tkinter as Tk
from pylab import *
from struct import *
from math import pow, log10
from sys import argv, exit

class srfclass:
    def __init__(self, fname, nx, nz, stkip, wstep, verbose=False):
        self.fname = fname
        self.nx = nx
        self.nz = nz
#self.dh = dh
#        self.dt = dt
#        self.nt = nt
        self.tskip = tskip
        self.wstep = wstep
        self.verbose = verbose

    def readsrf(self, ti):
        self.ti = ti
        askip = self.tskip * self.wstep
        self.V = zeros((self.nz, self.nx), dtype=float)
        
        if self.tskip == 0:
            cfname = self.fname
            fid = open(cfname, 'r')
        elif self.wstep == 0:
            cfname = self.fname
            fid = open(cfname, 'r')
            nskip = (ti // self.tskip) - 1
            if self.verbose:
                print("loading entry {0} in {1}".format(nskip, cfname))
            bskip = calcsize("<%df" % nskip) * self.nx * self.nz
            fid.seek(bskip, 0)
        
        fmt = "<%df" % self.nx
        sz = calcsize(fmt)

        for ii in range(self.nz):
            sr = fid.read(sz)
            self.V[ii,:] = unpack(fmt, sr)
        
        fid.close()

if len(argv) != 2:
    print("Usage: {0} directory of slip rate".format(argv[0]))
    exit()

dirname = argv[1]
nx = 2960
nz = 180
nt = 2000
dt = 0.06
wstep = 0
tskip = 1

srf = srfclass("%s" % dirname, nx, nz, tskip, wstep, verbose=True)

root = Tk.Tk()
root.title("%s" % dirname.split('/')[-1].upper())
root.protocol("WM_DELETE_WINDOW", root.quit)

tlabel = Tk.Label(root)
tlabel.pack()
ts = Tk.IntVar(root)
ts.set(1000)

fig = figure(figsize=(8,8))

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=Tk.TRUE)
toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.pack(side="top")

def replot(event):
    fig.clf()
    t = ts.get()
    fig.add_subplot(111)
    srf.readsrf(t)

    print("min, max slip rate: {0}, {1}".format(srf.V.min(), srf.V.max()))
    smax = abs(srf.V).max()
    imshow(flipud(srf.V), vmin=0, vmax=smax, aspect='auto')

    colorbar(orientation="horizontal")
    title("{variable} ({time} sec)".format(
                variable=dirname.split('/')[-1][-1].upper(), 
                time = "%5.3f" % (t * dt)))
    savefig("snap_{variable}_{time}s.pdf".format(
                      variable=dirname.split('/')[-1][-1].upper(),
                      time = "%5.3f" % (t * dt)))
#axis("image")
#    ax = gca()
#    ax.set_aspect(1)
    canvas.draw()

slider = Tk.Scale(root, from_=tskip, to=(nt - tskip), resolution=tskip,
                  orient=Tk.HORIZONTAL, variable=ts)
slider.pack(side="left", expand=Tk.TRUE, fill=Tk.X)
slider.bind('<ButtonRelease>', replot)
replot("<ButtonRelease>")

root.mainloop()

