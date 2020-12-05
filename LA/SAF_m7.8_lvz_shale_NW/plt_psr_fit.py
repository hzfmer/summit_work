#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 16:14:11 2018

@author: zhh076

To plot PSR-depth profile
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import os
from os.path import dirname as up

params = {'axes.labelsize': 24,
          'axes.titlesize': 24,
          'font.size': 24,
          'legend.fontsize': 20,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'xtick.major.size': 10,
          'xtick.minor.size': 5,
          'xtick.major.pad': 16,
          'ytick.major.size': 10,
          'ytick.minor.size': 5,
          'ytick.major.pad': 16}
rcParams.update(params)
np.seterr(divide='ignore', invalid='ignore')

nx = 2960
nz = 180
dx = 0.1

rock =os.getcwd().split('/')[-1][13:]
if "sandstone" in rock:
    r0_u = 0.4
    r0_w = 0.5
    GSI = 50
elif "shale" in rock:
    r0_u = 0.3
    r0_w = 0.4
    GSI = 30

u_loc = GSI * 3 / 50
w_loc = u_loc + 0.5

ofname = "PSR_depth_%s" % rock
ofile = open('%s.txt' % ofname,'w')

for ii in [chr(xx) for xx in range(ord("A"), ord("J"))]:
    try:
        ratu_lin = np.fromfile(os.path.join(os.getcwd().replace("shale", "sandstone"),
                              "lin_{0}/output_dyn/maxrate_x_01.0Hz.bin").format(ii), 
                              'f').reshape(nz, nx)
    except:
        continue

    ratw_lin = np.fromfile(os.path.join(os.getcwd().replace("shale", "sandstone"), 
                          "lin_{0}/output_dyn/maxrate_z_01.0Hz.bin").format(ii), 
                          'f').reshape(nz, nx)
    try:
        ratu_nl = np.fromfile("nl_{0}/output_dyn/maxrate_x_01.0Hz.bin".format(ii), 
                              'f').reshape(nz, nx)
    except:
        continue
    ratw_nl = np.fromfile("nl_{0}/output_dyn/maxrate_z_01.0Hz.bin".format(ii), 
                          'f').reshape(nz, nx)
    
    ratu_lin_dep = flipud(mean(ratu_lin, axis=1))
    ratw_lin_dep = flipud(mean(ratw_lin, axis=1))
    ratu_nl_dep = flipud(mean(ratu_nl, axis=1))
    ratw_nl_dep = flipud(mean(ratw_nl, axis=1))
    div_u = ratu_nl_dep / ratu_lin_dep
    div_u[div_u == 0] = 1
    div_u[isinf(div_u)] = 1
    div_w = ratw_nl_dep / ratw_lin_dep
    div_w[div_w == 0] = 1
    div_w[isinf(div_w)] = 1

    id_w = argwhere(absolute(diff(ratw_lin_dep)) < 0.03)
    id_w = id_w[2]
    id_u = 1 * id_w
    fid_stress = open(os.path.join(up(os.getcwd()),
                                        'source_7.8_so_cvm_lvz/source_%s.txt' % ii))
    stress = float(fid_stress.readlines()[-1].strip('\n').split(':')[-1])

    delta_r = (stress - 7e6) / 7e6 / 10
    r_u = 0.95
    r_w = 0.8
    r_u = r_u - delta_r
    r_w = r_w - delta_r
    ofile.write('{0}\nr_u: {1}\nr_w: {2}\n'.format(ii,r_u,r_w))
    ofile.write('delta_r: {0}\n'.format(delta_r))
    ratio_u = r0_u + np.random.normal(0, 0.02, nz) + (r_u - r0_u) * (1 -
            exp(-linspace(0, u_loc * nz / id_u, nz)))
    ratio_w = r0_w + np.random.normal(0, 0.02, nz) + (r_w - r0_w) * (1 -
                exp(-linspace(0, w_loc * nz / id_w, nz)))
    ax = figure(figsize=(7,8)) 
    plot(div_u, arange(nz) * dx, color='k', linewidth=2.5,
         label='data (r0 = {0})'.format('%.2f' % div_u[0]))
    plot(ratio_u, arange(nz) * dx, 'r:', linewidth=3,
         label='fit (r0 = {0})'.format('%.2f' % ratio_u[0]))
    xlabel('ratio')
    ylabel('Depth (km)')
    xlim(left=0)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    legend(loc=3)
    title('{0}_PSR_u NL/LIN'.format(ii))
    savefig('PSR_u_ratio_{0}_{1}.pdf'.format(rock, ii))
    close()
    
    ax = figure(figsize=(7,8)) 
    plot(div_w, arange(nz) * dx, color='k', linewidth=2.5,
         label='data (r0 = {0})'.format('%.2f' % div_w[0]))
    plot(ratio_w, arange(nz) * dx, 'r:', linewidth=3,
         label='fit (r0 = {0})'.format('%.2f' % ratio_w[0]))
    xlabel('ratio')
    ylabel('Depth (km)')
    xlim(left=0)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    ratio_w = ratio_w[::-1]
    legend(loc=3)
    title('{0}_PSR_w NL/LIN'.format(ii))
    savefig('PSR_w_ratio_{0}_{1}.pdf'.format(rock, ii))
    close()
    
# tar the figures
ofile.close()
script = """
tm=`date`
sed -i '1s#^#$tm#' {0}.txt;
tar -cvf {0}.tar PSR*.pdf {0}.txt; 
rm -f PSR*.{{txt,pdf}}
""".format(ofname)
os.system("bash -c '%s'" % script)
    
    
    
    
    
    
    
    
    

