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
    n = nx * ny
    skip = 10
    if M > 6:
        skip = 15
    try:
        data_gmpe = np.fromfile(f'sa_GMPE_M{M}_T{T}.bin', dtype='f').reshape((-1, 3))
    except FileNotFoundError:
        return 0
    dist, sa_median_gmpe, sa_sigma_gmpe = (data_gmpe[:, i] for i in range(3))
    print(f"GPME_M{M}_T{T} loaded!\n") 
    sa_cb = np.fromfile(f'{home_dir}gp_021219_Scott_{M}_noplas_GPU_' + freq + 
            f'/gmrotD50_{1 / T:05.2f}Hz.bin', dtype='float32')
    print("Spectrum of Cybershake read!\n")
    
    sa_cb = sa_cb.reshape((ny, nx)) / 10  # Convert m/s to g 
    sa_median_gmpe = sa_median_gmpe.reshape((ny, nx)) 
    fig, ax = plt.subplots(figsize=(8, 8))
    plt.rc('image', origin='lower')
    # vmax={6.35:0.06, 7.35: 0.2, 8.45: 1.2}
    # c = ax.imshow(sa_cb / sa_median_gmpe, extent=[0, nx * dh, 0, ny * dh], 
    #        norm=mpl.colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-0.5, vmax=vmax[M]), cmap='RdBu_r') 
    c = ax.imshow(sa_cb / sa_median_gmpe, extent=[0, nx * dh, 0, ny * dh], 
           norm=mpl.colors.LogNorm(vmin=0.2, vmax=5), cmap=discrete_cmap(8, 'RdBu_r'))
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    cbar = plt.colorbar(c, ax=ax, label=f'{T}s-SA (g) Ratio', orientation='horizontal',
            format='%d', spacing='proportional') 
    cbar.ax.xaxis.set_minor_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
    cbar.ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(subs=(0.2, 0.4, 0.6, 0.8)))
    # cbar.formatter.minor_thresholds = (1, 0.2)
    # cbar.set_ticks(np.arange(np.log10(c.colorbar.vmin), np.log10(c.colorbar.vmax), 0.2))
    # cbar.set_ticklabels(np.arange(np.log10(c.colorbar.vmin), np.log10(c.colorbar.vmax), 0.2))
    plt.tight_layout()
    fig.savefig(f"gp_021219_Scott_M{M}_T{T}_noplas_GPU_" + freq + "_imshow.png", dpi=400, bbox_inches='tight', pad_inches=0.05)

    sa_cb = np.ravel(sa_cb)
    sa_median_gmpe = np.ravel(sa_median_gmpe)
    idx = np.argsort(dist)
    dist = dist[idx]
    sa_median_gmpe = sa_median_gmpe[idx]
    sa_sigma_gmpe = sa_sigma_gmpe[idx]
    high = 10 ** (log10(sa_median_gmpe) + np.abs(log10(sa_sigma_gmpe)))
    low = 10 ** (log10(sa_median_gmpe) - np.abs(log10(sa_sigma_gmpe)))
    high2 = 10 ** (log10(sa_median_gmpe) + 2 * np.abs(log10(sa_sigma_gmpe)))
    low2 = 10 ** (log10(sa_median_gmpe) - 2 * np.abs(log10(sa_sigma_gmpe)))
    # x_bins = np.linspace(dist.min() + eps, dist.max(), int(np.sqrt(n / skip ** 2)))
    # y_bins = np.linspace(sa_cb.min() + eps, sa_cb.max(), int(np.sqrt(n)))
    x_bins = np.logspace(np.log10(dist.min() + eps), np.log10(dist.max()), 
            int(np.sqrt(n / skip ** 2)))
    y_bins = np.logspace(np.log10(sa_cb.min() + eps), np.log10(sa_cb.max()),
            int(np.sqrt(n / skip ** 2)))
    H, xedges, yedges = np.histogram2d(dist[::skip], sa_cb[::skip], bins=[x_bins, y_bins])
    fig, ax = plt.subplots(figsize=(8, 8))
    c = ax.pcolormesh(xedges, yedges, H.T, cmap=plt.get_cmap('twilight')) 
    c.set_rasterized(True)
    ax.scatter(dist[::skip], sa_median_gmpe[::skip], 2, 'r', 'o', rasterized=True)
    c1 = ax.fill_between(dist[::skip ** 2], low[::skip ** 2], high[::skip ** 2], 
            color='gray', alpha = 0.6, rasterized=True)
    # c2 = ax.fill_between(dist[::skip ** 2], low2[::skip ** 2], high2[::skip ** 2], 
    #        color='gray', alpha=0.3, rasterized=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("$R_{jb}$ distance (km)")
    ax.set_ylabel(f"SA {int(T)}s (g)")
    ax.set_xlim((1, dist.max()))
    plt.tight_layout()
    if M == 6.35:
        ax.set_ylim(bottom=1e-4)
    if M == 7.35:
        ax.set_ylim(bottom=5e-4)
    if M == 8.45:
        ax.set_ylim(bottom=5e-3)
    print(mpl.rcParams['agg.path.chunksize'])
    fig.savefig(f"gp_021219_Scott_M{M}_T{T}_noplas_GPU_" + freq + "_bin.png", dpi=400, bbox_inches='tight', pad_inches=0.05)


home_dir = "/ccs/home/hzfmer/scratch/LA/gp_rupture_test/LA/gp_rupture_test/"
nx = 6320
ny = 4200
dh = 0.1
eps = np.finfo(np.float32).eps
for M in [6.35]:#, 7.35, 8.45]:
    for T in [1, 2, 3, 5]:
        #plot_GMPE_CB(M, T, "1hz")
        # imshow_GMPE_CB(M, T, "1hz")
        plot_GMPE_CB(M, T, "2hz")
        # imshow_GMPE_CB(M, T, "2hz")
        plt.close('all')

