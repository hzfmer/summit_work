import numpy as np

nx, ny, nt = 540, 160, 800
vz = np.fromfile('output_sfc/SZ_0_0008000', dtype='float32').reshape(nt, ny, nx)
 
print(np.max(vz))
