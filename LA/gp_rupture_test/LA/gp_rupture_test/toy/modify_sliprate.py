import numpy as np

nx, nz, nt = 354, 212, 400
sr = np.fromfile('sliprate.bin_0', dtype='float32').reshape(-1, nz, nx)

sr = sr[:nt, : :]
sr.tofile('sliprate.bin')


