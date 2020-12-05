#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 25 23:57:50 2018

@author: zhh076
"""

from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt
from scipy.fftpack import fft
import os
from os.path import dirname as up
from filter_BU import filt_B

plt.rcParams['lines.linewidth'] = 1.5
def get_stat_name(path):
    f_list = os.listdir(path)
    dname_list = []
    for ii in f_list:
        if os.path.splitext(ii)[1] == '.bbp':
            dname_list.append(ii.strip('.bbp'))
    return dname_list


      
nstat = 351
tmax = 120
nt_lg = 12000
nt_rwg = 120000
dt_lg = tmax / nt_lg
dt_rwg = tmax / nt_rwg
dir_bbp = "/ccs/home/hzfmer/scratch/high_f/data/RWG_bbp/"

f_rwg = 1 / dt_rwg
f_lg = 1 / dt_lg 

fx = "/ccs/home/hzfmer/scratch/high_f/la_habra_large_gpu_dm_abc50/seis_x0120000"
fy = "/ccs/home/hzfmer/scratch/high_f/la_habra_large_gpu_dm_abc50/seis_y0120000"
fz = "/ccs/home/hzfmer/scratch/high_f/la_habra_large_gpu_dm_abc50/seis_z0120000"
sx = np.fromfile(fx,'f').reshape(nstat, nt_lg)
sy = np.fromfile(fy,'f').reshape(nstat, nt_lg)
sz = np.fromfile(fz,'f').reshape(nstat, nt_lg)
    
stat = np.zeros((nstat, 3, nt_rwg))
dname_list = get_stat_name(dir_bbp)
for i in range(nstat):
 #   dat[F.strip('.dat')] = fromfile(F,'b').reshape(len(dname_list),-1)  
    print(f"Reading station {i} / {nstat}\r", end="", flush=True)
    stat[i,:,:] = np.loadtxt(os.path.join(dir_bbp, '{0}.bbp'.format(dname_list[i]))
                        ,dtype='float').T

t_lg = np.linspace(0, tmax, nt_lg, endpoint=False)
t_rwg = np.linspace(0, tmax, nt_rwg, endpoint=False)

# for ii in range(nstat):
#     sname = dname_list[ii]
#     data = np.vstack((t_lg, sy[ii, :], sx[ii, :], sz[i, :])).astype('f')
#     np.savetxt(os.path.join(up(os.getcwd()), f'results/synthetics/{sname}.dat'), data.T)
    
for ii in range(nstat):
    sname = dname_list[ii]
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize = (10,8))
    plt.suptitle("Station: {0}".format(sname), fontsize=20)
    ax1.plot(t_rwg, stat[ii, 1, :], 'k')
    ax1.plot(t_lg, sx[ii,:], 'r', linestyle='--', linewidth=1)
    ax1.set_ylabel('V040', fontsize=16)
    ax1.legend(['RWG', 'AWP'],loc=4, fontsize=16)
    ax2.plot(t_rwg, stat[ii, 0, :], 'k')
    ax2.plot(t_lg, sy[ii,:], 'r', linestyle='--', linewidth=1)
    ax2.set_ylabel('V130 (m/s)', fontsize=16)
    ax3.plot(t_rwg, stat[ii, 2, :], 'k')
    ax3.plot(t_lg, sz[ii,:], 'r', linestyle='--', linewidth=1)
    ax3.set_ylabel('VZ', fontsize=16)
    ax3.set_xlabel('Time (s)', fontsize=16)
    fig.savefig(f'../results/{sname}_RWG_AWP_dm.png', dpi = 600,
                bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize = (10,8))
    plt.suptitle("Spectrum, Station: {0}".format(sname), fontsize=20)
    ax1.plot(np.linspace(0, f_rwg/2, int(nt_rwg/2)), 
             2 * dt_rwg * np.abs(fft(stat[ii, 1, :])[0 : int(nt_rwg/2)]), 'k')
    ax1.plot(np.linspace(0, f_lg/2, int(nt_lg/2)), 
             2 * dt_lg * np.abs(fft(sx[ii, :])[0 : int(nt_lg/2)]), 'r', linestyle='--', linewidth=1)
    ax1.set_ylabel('V040', fontsize=16)
    ax1.legend(['RWG', 'AWP'],loc=4, fontsize=16)
    
    ax2.plot(np.linspace(0, f_rwg/2, int(nt_rwg/2)), 
             2 * dt_rwg * np.abs(fft(stat[ii, 0, :])[0 : int(nt_rwg/2)]), 'k')
    ax2.plot(np.linspace(0, f_lg/2, int(nt_lg/2)), 
             2 * dt_lg * np.abs(fft(sy[ii, :])[0 : int(nt_lg/2)]), 'r', linestyle='--', linewidth=1)
    ax2.set_ylabel('V130', fontsize=16)
    
    ax3.plot(np.linspace(0, f_rwg/2, int(nt_rwg/2)), 
             2 * dt_rwg * np.abs(fft(stat[ii, 2, :])[0 : int(nt_rwg/2)]), 'k')
    ax3.plot(np.linspace(0, f_lg/2, int(nt_lg/2)), 
             2 * dt_lg * np.abs(fft(sz[ii, :])[0 : int(nt_lg/2)]), 'r', linestyle='--', linewidth=1)
    ax3.set_ylabel('VZ', fontsize=16)
    ax3.set_xlabel('Frequency (Hz)', fontsize=16)
    
    plt.xlim([0, 5])
    fig.savefig(f'../results/{sname}_RWG_AWP_spec_dm.png', dpi = 600,
                bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
        
    
