#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import struct
import collections

mx, my, mz = 5940, 4860, 108
# The sources are located in the second block
dh0, dh1 = 8, 24
nsrc = 15625
nt_src = 6600
gbuf = 1100

# For the second block in Z direction, there are (mz - 7) fine grids above, and (mz - 8)dh distance
# For the second block in Y direction, shift one more fine grid further, so yoff is dh smaller
seen = collections.defaultdict(list)
xoff, yoff, zoff = 0, -dh0, (mz - 8) * dh0
with open("momrate.dat", 'rb') as fidin, \
        open("source_1", 'wb') as fout:
    for i in range(nsrc):
        ix, iy, iz = struct.unpack('3I', fidin.read(12))
        ix = round(((ix - 1) * dh0 - xoff) / dh1) + 1
        iy = round(((iy - 1) * dh0 - yoff) / dh1) + 1
        iz = round(((iz - 1) * dh0 - zoff) / dh1) + 1
        seen[ix, iy, iz] += i,
        idx = np.array([ix, iy, iz], dtype='int32')
        mom = np.frombuffer(fidin.read(nt_src * 6 * 4), dtype='float32')
        idx.tofile(fout)
        mom.tofile(fout)
print(len(seen), nsrc)

with open("source_1", "r+b") as fid:
    for k, v in seen.items():
        tmp = np.zeros((6 * nt_src,), dtype='float32')
        for i in v:
            fid.seek(i * (6 * 4 * nt_src + 12), 0)
            ix, iy, iz = struct.unpack('3I', fid.read(12))
            assert (ix, iy, iz) == k
            tmp += np.frombuffer(fid.read(6 * 4 * nt_src), dtype='float32')
        for i in v:
            fid.seek(i * (6 * 4 * nt_src + 12) + 12, 0)
            tmp.tofile(fid)
try:
    open('source_0', 'ab', 0).close()
except OSError:
    pass
