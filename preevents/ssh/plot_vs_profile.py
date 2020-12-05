#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='paper', style='whitegrid', font_scale=1, rc={"xtick.bottom" : True, "ytick.left" : True})



plt.tight_layout(rect=[0, 0.03, 1, 0.95])
rcparams = {'font.size': 16,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 14,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'lines.linewidth': 2}
plt.rcParams.update(rcparams)

nz, ny, nx, nvar = 400, 400, 400, 5

mesh = ['mesh_TKCH05_081318.bin',
        'mesh_TKCH05_052119_2v680_bd98.bin']

labels = ['Borehole', 'Simplified']
colors = ['#ff7f0e', '#2ca02c']
fig, ax = plt.subplots(figsize=(6,6))
plt.gca().invert_yaxis()
fig2, ax2 = plt.subplots(figsize=(6,6))
plt.gca().invert_yaxis()
for i in range(len(mesh)):
    dat = np.fromfile(mesh[i], dtype='float32').reshape((nz,ny, nx, nvar))
    vs = dat[:, 201, 200, 1]
    #ax.plot(vs, np.arange(nz) * 2.5, c=colors[i], label=labels[i]) 
    #ax2.plot(vs[0:40], np.arange(40) * 2.5, c=colors[i], label=labels[i]) 
    ls = '--' if i == 1 else '-'
    ax.plot(vs, np.arange(nz) * 2.5, ls=ls, label=labels[i]) 
    ax2.plot(vs[0:40], np.arange(40) * 2.5, ls=ls, label=labels[i]) 

ax.set_xlabel('Vs (m/s)')
ax.set_ylabel('Depth (m)')
ax.legend()
fig.savefig(f"TKCH05_vs_profiles_052119.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

ax2.set_xlabel('Vs (m/s)')
ax2.set_ylabel('Depth (m)')
ax2.legend()
fig2.savefig(f"TKCH05_vs_profiles_shallow_052119.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
#v_pre = 0
#for i in range(nz):
#    if vs[i] != v_pre:
#        print(i, vs[i])
#        v_pre = vs[i]
#

def fetch_TKCH05(depth):
    if depth <= 6:
        return 140
    elif depth <= 18:
        return 430
    elif depth <= 30:
        return 660
    elif depth <= 80:
        return 770
    elif depth <= 98:
        return 640
    elif depth <= 295:
        return 1100
    
def fetch_HKD090(depth):
    if depth <= 2:
        return 137
    elif depth <= 5:
        return 97
    else:
        return 855
    
def plot_discret_curve(depths, vs=vs):
    
    v_HKD090 = list(map(fetch_HKD090, np.arange(0, 15)))
    fig, ax = plt.subplots(dpi=600)
    for depth in depths:
        v_TKCH05 = list(map(fetch_TKCH05, depth))
        inc = depth[1] - depth[0]
        ax.plot(v_TKCH05, depth, ls='-', label=f'TKCH05')
    ax.plot(v_HKD090, np.arange(0, 15), ls='-', label='HKD090')
    ax.plot(vs[:48], np.arange(48) * 2.5, ls='-.', label='Simplified')
    ax.set_xlabel('Vs (m/s)')
    ax.set_ylabel('Depth (m)')
    ax.legend(loc=3)
    ax.invert_yaxis()

vs[40] = 1100
plot_discret_curve([np.arange(0, 120, d) for d in [2.5]])
