#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
'''
Max and min of the mesh.
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import re

dir_name = os.getcwd().split('/')[-1]
if 'large' in dir_name:
    X = 180000
    Y = 135000
    Z = 61440
    dx = 20
elif 'small' in dir_name:
    X = 28000
    Y = 28000
    Z = 12000
    dx = 8
else:
    print("Please specify resolution in the directory name\n")
    sys.exit(-1)

dx = int(re.findall(r'\d+m', dir_name)[0][:-1])
    
nx = X // dx
ny = Y // dx
nz = Z // dx

fmax = 5
ppw = 8

try:
    vp = np.genfromtxt('vp_depth.dat')
    vpmax = vp[:, 0]
    vpmin = vp[:, 1]
    print(f"vpmax={np.max(vpmax[:108])}, {np.max(vpmax[:1250])}, vpmin={np.min(vpmin)}")
    vs = np.genfromtxt('vs_depth.dat')
    vsmax = vs[:, 0]
    vsmin = vs[:, 1]
    print(f"vsmax={np.max(vsmax[:108])}, {np.max(vsmax[:1250])}, vsmin={np.min(vsmin)}")
except IOError:
    with open('mesh', 'rb') as fidin:
        vsmax = np.zeros((nz, ))
        vsmin = np.zeros((nz, ))
        vpmax = np.zeros((nz, ))
        vpmin = np.zeros((nz, ))
        for ii in range(nz):
            buf = fidin.read(12 * nx * ny)
            dat = np.frombuffer(buf, dtype='f', count=3*nx*ny).reshape(-1,3)[:,0:2]
            vpmax[ii] = dat[:, 0].max()
            vpmin[ii] = dat[:, 0].min()
            vsmax[ii] = dat[:, 1].max()
            vsmin[ii] = dat[:, 1].min()
            print("\rReading layer {0}/{1}".format(ii,nz))
        vp=np.vstack((vpmax, vpmin)).T
        np.savetxt('vp_depth.dat',vp, fmt='%f')
        vs=np.vstack((vsmax, vsmin)).T
        np.savetxt('vs_depth.dat',vs, fmt='%f')

z0 = np.argwhere(vsmin > dx * fmax * ppw)[0][0] * dx / 1000
z1 = np.argwhere(vsmin > dx * fmax * 3 * ppw)[0][0] * dx / 1000
z2 = np.argwhere(vsmin > dx * fmax * 9 * ppw)[0][0] * dx / 1000
n1 = int(nx * ny * (z1 * 1000 // dx + 7))
n2 = int(nx * ny * ((z2 - z1) * 1000 // dx // 27 + (7 - 1) * 3 + 1))
n3 = int(nx * ny * (nz - z2 // dx) / 729)

print(f"Vmin = {min(vsmin)}\n")
print(f"Smallest ppw = {min(vsmin) / dx / fmax}\n")
print(f"First overlap at {z1 * 1000 // dx + 7}\n")  # Overlap
print(f"Second overlap at {z2 * 1000 // dx + (7 - 1) * 3 + 1}\n")  # Overlap
print(f'Total number of grids:\n'
      f'Layer 1, {n1}\n'
      f'Layer 2, {n2}\n' 
      f'Layer 3, {n3}\n'
      f'Total  , {n1 + n2 + n3} = {n1 + n2 + n3:e}\n')

fig, ax = plt.subplots()
ax.invert_yaxis()
ax.plot(vsmin/fmax/dx, np.arange(nz) / 1000 * dx, linewidth=2, label='vmin / dx / fmax')
ax.axhspan(0, z0, color='g', alpha=0.5, label='dx < 8')
ax.axhspan(z0, z1, color='b', alpha=0.5, label=f'Fine, z={z1}km')
ax.axhspan(z1, z2, color='r', alpha = 0.5, label=f'Medium, z={z2}km')
ax.axhspan(z2, nz * dx // 1000, alpha = 0.5, label='Coarse')
#ax.plot(np.full((nz,), ppw), np.arange(nz) * dx / 1000, '--', linewidth=3, label=f'ppw={ppw}')
#ax.plot(np.full((nz,), 3*ppw), np.arange(nz) * dx / 1000, '--', linewidth=3, label=f'ppw={3 * ppw}, depth = {z1 / 1000}km')
#ax.plot(np.full((nz,), 9*ppw), np.arange(nz) * dx / 1000, '--', linewidth=3, label=f'ppw={9 * ppw}, depth = {(z1 + z2) / 1000}km')
ax.legend()
ax.set_xlabel('PPW')
ax.set_ylabel('Depth (km)')
fig.savefig('ppw_depth.png', dpi=600, bbox_inches='tight', pad_inches=0.1)

