'''Compute suggested mesh split
Input:
    dh              Finest grid spacing
    [z1, z2, ...]   Expected depths of the bottom of overlapping zones

Reference:
    mesh_maxmin.py  Locate overlapping zones respecting ppw requirements
'''
import numpy as np
import sys

BLOCK_SIZE_Z = 4
z = []
dh = float(sys.argv[1])
for i in range(2, len(sys.argv)):
    z.append(float(sys.argv[i]))

h = {}
nz = {}
for i in range(len(z)):
    h[i] =  dh * 3 ** i
    if i < len(z) - 1:
        nz[i] = int(np.ceil(z[i] - sum((nz[j] - 8) * h[j] for j in range(i))) / h[i]) + 1 
        nz[i] = int(np.ceil(nz[i] / BLOCK_SIZE_Z)) * BLOCK_SIZE_Z
    else:
        nz[i] = int(np.floor(z[i] - sum((nz[j] - 8) * h[j] for j in range(i))) / h[i]) + 1
        nz[i] = nz[i] // BLOCK_SIZE_Z * BLOCK_SIZE_Z
    print(f"Depth preset at {z[i]}; real depth at "
          f"{sum((nz[j] - 1 - 7 * (j < i))* h[j] for j in range(i + 1))}\n")

print(f"Set the number of layers to be:\n"
      f"{list(nz.values())}")

