import numpy as np
from scipy.signal import resample


nsrc = 15625
nt = 12500
ncomp = 6
nbit = 4
dh_0, dh_1 = 8, 24
dtop = (160 - 8) * dh_0
x_off, y_off = -1000, -1000
with open("source_1", "wb") as fid, \
    open("../dhyp1.50_s372823598_q100f00_orig_vs200/source_1", "rb") as fid_in:
        for _ in range(nsrc):
            ix, iy, iz = np.fromfile(fid_in, dtype='int32', count=3)
            # ix = int((ix - 1) * dh_1 / dh_0) + x_off
            # iy = int((iy - 1) * dh_1 / dh_0) + y_off
            # iz = int(((iz - 1) * dh_1 + dtop) / dh_0)
            ix = ix // 5
            iy = iy // 5
            iz = iz 
            idx = np.array([ix, iy, iz], dtype='int32')
            idx.tofile(fid)
            data = np.fromfile(fid_in, dtype='float32', count=nt * ncomp).reshape(nt, ncomp)
            data.tofile(fid)

