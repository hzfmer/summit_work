#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

nx, ny = 6320, 4200
M = 6.35
T = 1

data_gmpe = np.fromfile(f'sa_GMPE_M{M}_T{T}.bin', dtype='f').reshape((ny, nx, 3))
dist2 = data_gmpe[:, :, 0]
dist = np.fromfile(f'mesh_rjb_M{M}.bin', dtype='float32').reshape((ny, nx))
fig, ax = plt.subplots()
ax.imshow(dist-dist2)
fig.savefig(f"delta_dist", dpi=600, bbox_inches='tight', pad_inches=0.05)
print(np.max(dist-dist2))
