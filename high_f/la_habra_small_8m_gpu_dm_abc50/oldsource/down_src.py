#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import struct

nx, ny = 3500, 3500
mx, my, mz = nx // 54 * 54, ny // 54 * 54, 108
px, py = (nx - mx) // 2, (ny - my) // 2
# The sources are located in the second block
dh, dh0, dh1 = 8, 20, 24
nsrc = 1560
nt = 5000
gbuf = 1000
xoff, yoff, zoff = px * dh, py * dh, 0

with open("momrate_LaHabra.zf.100m.rev.bin", 'rb') as fidin, \
        open("source_0", 'wb') as fout:
    for i in range(nsrc):
        ix, iy, iz = struct.unpack('3I', fidin.read(12))
        ix = round(((ix - 1) * dh0 - xoff) / dh) + 1
        iy = round(((iy - 1) * dh0 - yoff) / dh) + 1
        iz = round(((iz - 1) * dh0 - zoff) / dh) + 1
        idx = np.array([ix, iy, iz], dtype='int32')
        mom = np.frombuffer(fidin.read(nt * 6 * 4), dtype='float32')
        idx.tofile(fout)
        mom.tofile(fout)

# For the second block in Z direction, there are (mz - 7) fine grids above, and (mz - 8)dh distance
# For the second block in Y direction, shift one more fine grid further, so yoff is dh smaller
xoff, yoff, zoff = px * dh, (py - 1) * dh, (mz - 8) * dh
with open("momrate_LaHabra.zf.100m.rev.bin", 'rb') as fidin, \
        open("source_1", 'wb') as fout:
    for i in range(nsrc):
        ix, iy, iz = struct.unpack('3I', fidin.read(12))
        ix = round(((ix - 1) * dh0 - xoff) / dh1) + 1
        iy = round(((iy - 1) * dh0 - yoff) / dh1) + 1
        iz = round(((iz - 1) * dh0 - zoff) / dh1) + 1
        idx = np.array([ix, iy, iz], dtype='int32')
        mom = np.frombuffer(fidin.read(nt * 6 * 4), dtype='float32')
        idx.tofile(fout)
        mom.tofile(fout)
