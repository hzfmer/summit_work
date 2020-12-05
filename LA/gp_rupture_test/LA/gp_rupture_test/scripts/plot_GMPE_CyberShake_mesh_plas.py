#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue March  19 18:37:05 2019

@author: zhh076
"""

import numpy as np
from numpy import log10, exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
sys.path.append("/ccs/home/hzfmer/scratch/LA/gp_rupture_test/LA/gp_rupture_test/scripts")
import BSSA_2014_nga
#matplotlib.rcParams.update({'font.size': 20})

import matplotlib.style as mplstyle
mplstyle.use('fast')
mpl.rc('font', size=16)

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

           
def plot_GMPE_CB(M, T, freq):
    nx, ny = 6320, 4200
    sa_tmp = np.fromfile(f'{home_dir}gp_021219_Scott_{M}_plas_' + freq + 
            f'/gmrotD50_{1 / T:05.2f}Hz.bin', dtype='float32')
    sa_cb = sa_tmp / 10
    data_gmpe = np.fromfile(f'sa_GMPE_M{M}_T{T}.bin', dtype='f').reshape((-1, 3))
    dist, sa_median_gmpe, sa_sigma_gmpe = (data_gmpe[:, i] for i in range(3))
    
    skip = 10
    dist, sa_cb, sa_median_gmpe, sa_sigma_gmpe = (x.reshape((ny, nx)) for x in \
            (dist, sa_cb, sa_median_gmpe, sa_sigma_gmpe))
    dist, sa_cb, sa_median_gmpe, sa_sigma_gmpe = (np.ravel(x[::skip, ::skip]) \
            for x in (dist, sa_cb, sa_median_gmpe, sa_sigma_gmpe))
    #idx = np.argsort(dist)
    #dist = dist[idx]
    #sa_cb = sa_cb[idx]
    #sa_median_gmpe = sa_median_gmpe[idx]
    #sa_sigma_gmpe = sa_sigma_gmpe[idx]
    high = 10 ** (log10(sa_median_gmpe) + np.abs(log10(sa_sigma_gmpe)))
    low = 10 ** (log10(sa_median_gmpe) - np.abs(log10(sa_sigma_gmpe)))

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.errorbar(dist, sa_median_gmpe, yerr=high - sa_median_gmpe, marker='o', 
            lolims=True, color='lightskyblue', linestyle='None', errorevery=5*skip,
            label='BSSA14')
    ax.scatter(dist[::5*skip], sa_cb[::5*skip], marker='D', color='orangered', label='AWP')
    # ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("$R_{jb}$ distance (km)")
    ax.set_ylabel(f"SA-{int(T)}s (g)")
    ax.set_xlim((1, dist.max()))
    plt.legend(loc=1)
    plt.tight_layout()
    if M == 6.35:
        ax.set_ylim(bottom=1e-4)
    if M == 7.35:
        ax.set_ylim(bottom=5e-4)
    if M == 8.45:
        ax.set_ylim(bottom=5e-3)
    fig.savefig(f"gp_021219_Scott_M{M}_T{T}_plas_" + freq + "_mesh_vs30.png", dpi=600, bbox_inches='tight', pad_inches=0.05)


home_dir = "/ccs/home/hzfmer/scratch/LA/gp_rupture_test/LA/gp_rupture_test/"
nx = 6320
ny = 4200
dh = 0.1
eps = np.finfo(np.float32).eps
for M in [6.35, 7.35, 8.45]:
    for T in [1, 2, 3, 5]:
        plot_GMPE_CB(M, T, "1hz")
        plot_GMPE_CB(M, T, "2hz")
        plt.close('all')


