#! /usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
from math import pow, log10
from matplotlib.pylab import *
import os

params={'axes.labelsize': 24,
        'axes.titlesize': 24,
        'font.size': 24,
        'legend.fontsize': 20,
        'xtick.labelsize': 24,
        'ytick.labelsize': 24,
        'xtick.major.size':10,
        'xtick.minor.size':5,
        'ytick.major.size':10,
        'ytick.minor.size':5}
rcParams.update(params)
#rc('font', family='serif')
#rc('font', serif='Times New Roman')
def mom2mag(m0):
    mw = 2. / 3. * (log10(m0) - 9.1)
    return mw

nx = 2960
nxx = 3000
nz = 180
xoff = 20
dx = 100

for vdir in [x + y + '/' for x in ['nl_', 'lin_'] for y in [chr(xx) for xx in range(ord("A"), 1+ord("D"))] + ['lin/', 'nl/']]:
    try:
        slipu = np.fromfile(vdir + "output_dyn/slipu", 'f', count = nx * nz).reshape(nz, nx)
    except Exception as e:
        continue

    slipw = np.fromfile(vdir + "output_dyn/slipw", 'f', count = nx * nz).reshape(nz, nx)
    rtime = np.fromfile(vdir + "output_dyn/rtime", 'f', count = nx * nz).reshape(nz, nx)
    ratu = np.fromfile(vdir + "output_dyn/maxrate_x_01.0Hz.bin", 'f', count = nx * nz).reshape(nz, nx) 
    ratw = np.fromfile(vdir + "output_dyn/maxrate_z_01.0Hz.bin", 'f', count = nx * nz).reshape(nz, nx) 
    rat = np.sqrt(ratu ** 2 + ratw ** 2)
    vp, vs, rho = np.loadtxt("../source_7.8_so_cvm_lvz/fault_3000.dat", unpack=True, usecols=[3,4,5])
    vp = vp.reshape(nz, nxx)[:, xoff : nxx - xoff]
    vs = vs.reshape(nz, nxx)[:, xoff : nxx - xoff]
    rho = rho.reshape(nz, nxx)[:, xoff : nxx - xoff]
    mu = rho * vs * vs
    
    slip = np.zeros((nz, nx), dtype = float)
    m0 = 0.
    for k in range(nz):
        for l in range(nx):
            slip[k, l] = np.sqrt(pow(slipu[k,l], 2.) + pow(slipw[k, l], 2.))
            m0 += mu[k, l] * pow(dx, 2.) * slip[k, l]
    
    mag = mom2mag(m0)
    
    x = np.linspace(0, (nx - 1) * dx, nx) / 1.e3
    z = np.flipud(np.linspace(0, (nz - 1) * dx, nz)) / 1.e3
    
    smax = np.ceil(slip.max())
    fig = figure(figsize=(12.6, 4))
    ax = subplot(111)
    cnt = np.linspace(0, smax, 21)
#print(cnt)
    contourf(x, z, abs(slipu), cnt, cmap=cm.jet)
    cb = colorbar(orientation = "vertical", shrink = 0.7, ticks = np.arange(0., smax + 0.1, 2.), aspect = 8)
    cb.set_label("Slip (m)")
    ax.get_yaxis().set_ticks(np.arange(0., 18., 5.))
    
    cnt = np.arange(0., 100., 1.)
    ct = contour(x, z, rtime, cnt, colors='k')
    axis((0., 296., 18., 0.))
    xlabel("Strike (km)")
    ylabel("Dip (km)")
    
#    vdir = os.path.basename(os.getcwd())[-1]
    rock =os.getcwd().split('/')[-1][13:]
#title('{0}_{1} ($M_w = {2}$)'.format(rock, vdir.strip('/'), '%.2f' % mag))
    fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 1, top = 0.87)
#   savefig(vdir + 'slip_%s.png' % vdir.strip('/'))
    savefig('slip_{0}_{1}.png'.format(rock, vdir.strip('/')))
    plt.close('all')
    
#smax = np.ceil(ratu[100,:].max())
    smax = 8
    fig = figure(figsize=(12.6, 4))
#ax = subplot(111)
    cnt = np.linspace(0, smax, 21)
#print(cnt)
    axMap = plt.axes([0.1, 0.25, 0.85, 0.3])
    contourf(x, z, abs(rat), cnt, cmap=cm.jet)
#   cb = colorbar(orientation = "vertical", shrink = 0.7, ticks = np.arange(0., smax + 0.1, 2.), aspect = 8)
#cb.set_label("Slip rate(m/s)")
    axMap.get_yaxis().set_ticks(np.arange(0., 18., 5.))
    axMap.set_xlabel('Strike (km)')
    axMap.set_ylabel('Depth (km)')
    cnt = np.arange(0., 100., 1.)
    ct = contour(x, z, rtime, cnt, colors='k')
    axis((0., 296., 18., 0.))
    axTop = plt.axes([0.1, 0.55, 0.85, 0.25 ])
    axTop.plot(x, rat[-1, :])
    axTop.set_yticks(np.round([tmp * np.max(rat[-1, :]) for tmp in [0.5, 1.05]]))
    axTop.set_xlim(0, x[-1])
    axTop.yaxis.tick_right()
    axTop.set_xticks([])
#    vdir = os.path.basename(os.getcwd())[-1]
    rock =os.getcwd().split('/')[-1][13:]
#title('{0}_{1} ($M_w = {2}$)'.format(rock, vdir.strip('/'), '%.2f' % mag))
    fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 1, top = 0.87)
#   savefig(vdir + 'slip_%s.png' % vdir.strip('/'))
    savefig('ratu_{0}_{1}.png'.format(rock, vdir.strip('/')),dpi=600)
    plt.close('all')

#smax = np.ceil(ratw[10,:].max())
    smax = 6
    fig = figure(figsize=(12.6, 4))
    ax = subplot(111)
    cnt = np.linspace(0, smax, 21)
#print(cnt)
    contourf(x, z, abs(ratw), cnt, cmap=cm.jet)
    cb = colorbar(orientation = "vertical", shrink = 0.7, ticks = np.arange(0., smax + 0.1, 2.), aspect = 8)
    cb.set_label("Slip rate(m/s)")
    ax.get_yaxis().set_ticks(np.arange(0., 18., 5.))
    
    cnt = np.arange(0., 100., 1.)
    ct = contour(x, z, rtime, cnt, colors='k')
    axis((0., 296., 18., 0.))
    xlabel("Strike (km)")
    ylabel("Dip (km)")
    
#    vdir = os.path.basename(os.getcwd())[-1]
    rock =os.getcwd().split('/')[-1][13:]
#title('{0}_{1} ($M_w = {2}$)'.format(rock, vdir.strip('/'), '%.2f' % mag))
    fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 1, top = 0.87)
#   savefig(vdir + 'slip_%s.png' % vdir.strip('/'))
    savefig('ratw_{0}_{1}.png'.format(rock, vdir.strip('/')),dpi=600)
    plt.close('all')

script = """
tar -cvf ratu_{0}.tar ratu*.png; 
tar -cvf ratw_{0}.tar ratw*.png; 
tar -cvf slip_{0}.tar slip*.png; 
rm -f *.png
""".format(rock)
os.system("bash -c '%s'" % script)
