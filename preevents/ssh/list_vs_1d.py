#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

nz, ny, nx, nvar = 400, 400, 400, 5

mesh = ['mesh_TKCH05_092718_left.bin',
        'mesh_TKCH05_042919_A78_bd975.bin',
        'mesh_TKCH05_042919_A78_bd975_left.bin',
        'mesh_TKCH05_043019_A78_bd975.bin',
        'mesh_TKCH05_042919_A88_bd975.bin']
#mesh = ['mesh_TKCH05_het_042919_A78_bd975_h005l100.bin',
#        'mesh_TKCH05_het_042919_A78_bd975_left_h005l100.bin',
#        'mesh_TKCH05_het_043019_A78_bd975_h005l100.bin',
#        'mesh_TKCH05_het_A88_bd975_h005l100.bin']

fig, ax = plt.subplots(figsize=(6,6))
plt.gca().invert_yaxis()
fig2, ax2 = plt.subplots(figsize=(6,6))
plt.gca().invert_yaxis()
for i in range(len(mesh)):
    dat = np.fromfile(mesh[i], dtype='float32').reshape((nz,ny, nx, nvar))
    vs = dat[:, 201, 200, 1]
    ax.plot(vs + 5 * i, np.arange(nz) * 2.5, label=f'{i}') 
    ax2.plot(vs[0:40], np.arange(40) * 2.5, label=f'{i}') 

ax.set_xlabel('Vs (m/s)')
ax.set_ylabel('Depth (m)')
ax.legend()
fig.savefig(f"TKCH05_vs_profiles_043019.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

ax2.set_xlabel('Vs (m/s)')
ax2.set_ylabel('Depth (m)')
ax2.legend()
fig2.savefig(f"TKCH05_vs_profiles_shallow_043019.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
#v_pre = 0
#for i in range(nz):
#    if vs[i] != v_pre:
#        print(i, vs[i])
#        v_pre = vs[i]
#
