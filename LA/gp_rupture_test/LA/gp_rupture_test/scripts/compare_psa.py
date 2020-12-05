#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os

rcparams = {'font.size': 16,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 14,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'lines.linewidth': 1}
plt.rcParams.update(rcparams)


for M in [6.35, 7.35, 8.45]:
    fig, ax = plt.subplots(3, 1, figsize=(6,6))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.suptitle(f"M = {M}, gmrot / cb")
    data_gmpe = np.fromfile(f'/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/gp_rupture_test/scripts/sa_GMPE_site_M{M}_T1.bin', dtype='f').reshape((-1, 3))
    dist = data_gmpe[:, 0]
    ll = np.genfromtxt(f'stat_M{M}.txt', dtype='int', skip_header=1)
    nsite = len(ll)
    for i, T in enumerate([1, 2, 5]):
        sa_cb = np.fromfile(f'/ccs/home/hzfmer/scratch/LA/gp_rupture_test/LA/data_for_Zhifeng/sa_cb_site_M{M}_T{T}.bin', dtype='f')
        sa_tmp = np.fromfile(f'/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/gp_rupture_test/gp_021219_Scott_{M}_noplas_GPU_2hz/gmrotD50_{1/T:05.2f}Hz.bin', dtype='f').reshape((4200, 6320))
        sa_awp = [sa_tmp[ll[i][1], ll[i][0]] / 9.8 for i in range(nsite)]
        ratio = np.divide(sa_awp, sa_cb)
        ax[i].scatter(dist, ratio)
        ax[i].set_ylabel('ratio')
        ax[i].legend([f'T={T}s, avg={np.mean(ratio):.2f}'])

    ax[-1].set_xlabel('Rjb (km)')
    fig.savefig(os.getcwd() + "/" + f"psa_ratio_M{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
