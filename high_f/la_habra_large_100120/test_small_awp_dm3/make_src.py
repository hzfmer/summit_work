import numpy as np
from scipy.signal import resample


nsrc = 15625
nt = 10000
ncomp = 6
nbit = 4
dh_0, dh_1 = 8, 24
dtop = (160 - 8) * dh_0
x_off, y_off = -400, -1000
with open("source_1", "wb") as fid, \
    open("../test_run_dm/source_1", "rb") as fid_in:
        for _ in range(nsrc):
            ix, iy, iz = np.fromfile(fid_in, dtype='int32', count=3)
            # ix = int((ix - 1) * dh_1 / dh_0) + x_off
            # iy = int((iy - 1) * dh_1 / dh_0) + y_off
            # iz = int(((iz - 1) * dh_1 + dtop) / dh_0)
            ix = ix + x_off // 3
            iy = iy + y_off // 3
            iz = iz 
            idx = np.array([ix, iy, iz], dtype='int32')
            idx.tofile(fid)
            data = np.fromfile(fid_in, dtype='float32', count=nt * ncomp).reshape(nt, ncomp)
            data.tofile(fid)

