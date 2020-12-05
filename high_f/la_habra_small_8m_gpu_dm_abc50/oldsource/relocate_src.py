#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import struct
from scipy.signal import resample

nx, ny = 3500, 3500
mx, my = nx // 54 * 54, ny // 54 * 54
px, py = (nx - mx) // 2, (ny - my) // 2
dh = 20
nsrc = 1560
nt = 5000
gbuf = 1000
xoff, yoff, zoff = px + 1, py + 1, 1

dt0 = 0.001
dt1 = 0.001  # whether need resample, if dt changes
nt1 = int(nt * dt0 / dt1)
gbuf = int(gbuf * dt1 / dt0)

dir_src = "input/source.txt"
dir_src_xx = "input/source_xx"
dir_src_yy = "input/source_yy"
dir_src_zz = "input/source_zz"
dir_src_xz = "input/source_xz"
dir_src_yz = "input/source_yz"
dir_src_xy = "input/source_xy"

fxx = open(dir_src_xx, 'w')
fyy = open(dir_src_yy, 'w')
fzz = open(dir_src_zz, 'w')
fxz = open(dir_src_xz, 'w')
fyz = open(dir_src_yz, 'w')
fxy = open(dir_src_xy, 'w')

with open("momrate_LaHabra.zf.100m.rev.bin", 'rb') as fidin, open(f"{dir_src}", "w") as f_coord:
    f_coord.write(
            f"2.0.0\n"
            f"file=input/source\n"
            f"degree=3\n"
            f"steps={nt1}\n"
            f"stride=1\n"
            f"gpu_buffer_size={gbuf}\n"
            f"cpu_buffer_size={nt1 // gbuf}\n"
            f"length={nsrc}\n\n"
            f"coordinates\n")
    for i in range(nsrc):
        ix, iy, iz = struct.unpack('3I', fidin.read(12))
        mom = np.frombuffer(fidin.read(nt * 6 * 4), dtype='float32').reshape(nt, 6)
        f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n") 
        mom[:, 0].tofile(fxx)
        mom[:, 1].tofile(fyy)
        mom[:, 2].tofile(fzz)
        mom[:, 3].tofile(fxz)
        mom[:, 4].tofile(fyz)
        mom[:, 5].tofile(fxy)

        # resample(mom[:, 0], nt1).astype('float32').tofile(fxx)
        # resample(mom[:, 1], nt1).astype('float32').tofile(fyy)
        # resample(mom[:, 2], nt1).astype('float32').tofile(fzz)
        # resample(mom[:, 3], nt1).astype('float32').tofile(fxz)
        # resample(mom[:, 4], nt1).astype('float32').tofile(fyz)
        # resample(mom[:, 5], nt1).astype('float32').tofile(fxy)

fxx.close()
fyy.close()
fzz.close()
fxz.close()
fyz.close()
fxy.close()
print("Done.\n")



