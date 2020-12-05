#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import struct
from scipy.signal import resample
from post_processing.la_habra import read_param
import sys


# Read simulation parameters
tmax, dt, tskip, wstep, nfile = read_param('.') 

# Mesh info
mx, my = 3456, 3456
dh = 8
xoff, yoff, zoff = 0, 0, 0

x0, y0 = 1000, 800
width = 20
length = 20
depth = 10
nsite = 2 * depth * (width + length - 2) + width * length

# Hardcoded parameters
version = '2.0.0'
file_prefix = 'output_vol/receivers'
nsite = nsite  # Number of sites
steps = int(tmax / dt)
stride = tskip  # tiem steps to skip when writing
degree = 3  # interpolation
gpu_buffer_size = 500
cpu_buffer_size = 3
num_writes = 1 

dir_recv = "input/receiver.txt"


with  open(f"{dir_recv}", "w") as f_coord:
    f_coord.write(
            f"{version}\n"
            f"file={file_prefix}\n"
            f"length={nsite}\n"
            f"steps={steps}\n"
            f"stride={stride}\n"
            f"cpu_buffer_size={cpu_buffer_size}\n"
            f"gpu_buffer_size={gpu_buffer_size}\n"
            f"num_writes={num_writes}\n"
            f"degree={degree}\n\n"
            f"coordinates\n")

    for iz in range(depth + 1):
        if iz == depth:
            for j in range(width):
                for i in range(length):
                    ix, iy = x0 + i, y0 + j
                    f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n")
            break
        for i in range(length - 1):  # left bottom --> right bottom
            ix, iy = x0 + i, y0 + 0
            f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n")
        for j in range(width - 1):  # right bottom --> right top
            ix, iy = x0 + length - 1, y0 + j
            f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n")
        for i in range(length - 1, 0, -1):  # right top --> left top
            ix, iy = x0 + i, y0 + width - 1
            f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n")
        for j in range(width - 1, 0, -1):  # left top --> left bottom
            ix, iy = x0, y0 + j
            f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n")


assert steps % (cpu_buffer_size * gpu_buffer_size) == 0
print("Done.\n")



