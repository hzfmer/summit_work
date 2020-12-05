#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import struct

nx, ny = 9000, 6750
mx, my = 3500, 3500
px, py = (nx - mx) // 2, (ny - my) // 2
# The sources are located in the second block
nsrc = 1560
nt = 5000
xoff, yoff = px, py 

with open("momrate_LaHabra.zf.100m.rev.bin", 'rb') as fidin, \
        open("momrate_LaHabra.zf.100m.rev.large.bin", 'wb') as fout:
    for i in range(nsrc):
        ix, iy, iz = struct.unpack('3I', fidin.read(12))
        ix, iy , iz = ix + xoff, iy + yoff, iz
        idx = np.array([ix, iy, iz], dtype='int32')
        mom = np.frombuffer(fidin.read(nt * 6 * 4), dtype='float32')
        idx.tofile(fout)
        mom.tofile(fout)
