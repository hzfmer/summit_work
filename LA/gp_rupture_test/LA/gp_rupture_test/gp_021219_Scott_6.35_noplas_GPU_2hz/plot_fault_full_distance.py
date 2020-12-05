#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

nx, ny = 6320, 4200
dh = 0.1
for M in [6.35]:
    try:
        sites = np.genfromtxt('stat_M{M}.txt', skip_header=1)
    except OSError:
        sites = np.genfromtxt('stat.txt', skip_header=1)
    rjb_sites = np.fromfile(f'site_rjb_M{M}.bin', dtype='f')
    dist = np.fromfile(f'mesh_rjb_M{M}.bin', dtype='float32').reshape((ny, nx))
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(dist, extent=[0, nx * dh, 0, ny * dh], cmap='RdBu_r', origin='lower')
    for _ in range(5):
        i = np.random.randint(len(sites))
        ax.scatter(sites[i,0] * dh, sites[i, 1] * dh, label=f'{rjb_sites[i]:.2f}km')
    ax.legend()
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    cb = plt.colorbar(im, ax=ax, orientation='horizontal', shrink=0.6)
    cb.set_label('Rjb (km)')
    fig.savefig(f"mesh_rjb_M{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

