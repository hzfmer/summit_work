#!/sw/titan/python/anaconda3-5.1.0/bin/python
import numpy as np

nx = 1000
ny = 480
nz = 7
nt = 500
sx = 600
sy = 242
data_dir = "output_sfc/"
for F in ['X', 'Y', 'Z']:
    v = np.array([])
    for i in range(29):
        dat = np.fromfile(data_dir + f"S{F}_0_{(i + 1)*12500:07d}", 'f').reshape((nt, nz, ny, nx))
        v = np.append(v, dat[:,:,sy,sx])
    print(v.shape)
    v = v.reshape(-1, nz)
    np.savetxt(f"V{F}.txt", v)
