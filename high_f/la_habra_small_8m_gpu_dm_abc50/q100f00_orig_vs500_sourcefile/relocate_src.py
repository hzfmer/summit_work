#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import struct
from scipy.signal import resample
from post_processing.la_habra import read_param
import sys

def resample_src(stf, n, ntotal):
    """Resmaple src to length n, then append zeros in the end to length of ntotal
    
    """
    return np.pad(resample(stf, n), (0, ntotal - n)).astype('float32')

# Raw source parameters
nsrc = 15625
nt_raw = 5000
dt_raw = 0.001
tsrc = nt_raw * dt_raw  # length of source 
gbuf = 1000

# whether need resample, if dt_raw changes to dt_new
tmax, dt_new, tskip, _, _ = read_param('.') 
nt_new = int(nt_raw * dt_raw / dt_new)
ntotal = int(tmax / dt_new)

# Append zeros
print(f"TMAX = {tmax}, tsrc={tsrc}, nt_raw={nt_raw}, nt_new={nt_new}")
print(f"Upsample from {dt_raw} to {dt_new}, in total {ntotal} time steps")

# Mesh info
mx, my = 3456, 3456
dh = 8
xoff, yoff, zoff = 1, 1, 1

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

#with open("momrate.dat", 'rb') as fidin, \
#     open(f'fault_idx.txt', 'w') as f_idx:
#        for i in range(nsrc):
#            ix, iy, iz = struct.unpack('3I', fidin.read(12))
#            f_idx.write(f'{ix} {iy} {iz}\n')
#            fidin.seek(nt_raw * 6 * 4, 1)
#
#sys.exit(-1)

with open("../momrate.dat", 'rb') as fidin, \
     open(f'fault_idx.txt', 'w') as f_idx, \
     open(f"{dir_src}", "w") as f_coord:
    f_coord.write(
            f"2.0.0\n"
            f"file=input/source\n"
            f"degree=3\n"
            f"steps={ntotal}\n"
            f"stride=1\n"
            f"gpu_buffer_size={gbuf}\n"
            f"cpu_buffer_size={ntotal // gbuf}\n"
            f"length={nsrc}\n\n"
            f"coordinates\n")
    for i in range(nsrc):
        ix, iy, iz = struct.unpack('3I', fidin.read(12))
        f_idx.write(f'{ix} {iy} {iz}\n')
        mom = np.frombuffer(fidin.read(nt_raw * 6 * 4), dtype='float32').reshape(nt_raw, 6)
        f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n") 
        for j, f in enumerate([fxx, fyy, fzz, fxz, fyz, fxy]):
            resample_src(mom[:, j], nt_new, ntotal).tofile(f)

for f in [fxx, fyy, fzz, fxz, fyz, fxy]:
    f.close()
print("Done.\n")



