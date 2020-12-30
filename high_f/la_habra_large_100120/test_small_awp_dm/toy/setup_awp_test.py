#!/usr/bin/env python3

"""
A aggregrative script to set up a simple AWP test case with configurations.

Include topography. 
Use awp-TOPO's source/receiver implementation.

Feel free to adjust almost all manually assigned values.

Note 
----
Buffer:
    When a large block of data are read or written to disk, buffer is used for the memory sake.
    Take reading as example, each time CPU sends cpu_buffer_size of gpu_buffer_size to GPU, GPU reads gpu_buffer_size data each time, then flush, then next until all "cpu_buffer_size" data processed. Then CPU goes to next chunk of data.

Parameters:
steps/stride/length/num_writes:
`length` is the number of sources/receivers
`steps` is the total time step of source time function or simulation duration
`stride` is how many time steps are skipped when reading/writing. For source, it's almost always 1.
`num_writes` is the number of the chunks to save

Author
------
    Zhifeng Hu: SDSU/UCSD
    Email: hzfmer94@gmail.com
    04/29/20

"""


# Toy model

from pathlib import Path

import numpy as np


def generate_gp(tau, t):
    """

    A kinematic source generator GP2010, used in CyberShake.

    Input
    -----
    tau : float
        the duration of stf
    t : array of float
        time steps

    Note
    ----
    Do not modify any hard coded parameters here.

    """
    stf = np.zeros((len(t), 1))
    tau1 = 0.13 * tau
    tau2 = tau - tau1
    cn = np.pi / (1.5 * np.pi * tau1 + 1.2 * tau1 + 0.2 * np.pi * tau2)
    for i, tt in enumerate(t):
        if tt < 0:
            stf[i] = 0
        elif tt < tau1:
            stf[i] = cn * (0.7 - 0.7 * np.cos(np.pi * tt / tau1) +
                          + 0.6 * np.sin(0.5 * np.pi * tt / tau1))
        elif tt < 2 * tau1:
            stf[i] = cn * (1.0 - 0.8 * np.cos(np.pi * tt / tau1)
                          + 0.2 * np.cos(np.pi * (tt - tau1) / tau2))
        elif tt < tau:
            stf[i] = cn * (0.2 + 0.2 * np.cos(np.pi * (tt - tau1) / tau2))
        else:
            stf[i] = 0
    return stf


def resample_src(stf, n, ntotal):
    """Resmaple src to length n, then append zeros in the end up to a length of ntotal 
    """
    return np.pad(resample(stf.squeeze(), n), (0, ntotal - n)).astype('float32')


def gen_2d_gauss(nx, ny, mu_x=None, mu_y=None, sigma_x=None, sigma_y=None):
    """2D Gaussian topography
    Input
    -----
    (nx, ny) : (int, int)
        Size of the topography mesh
    (mu_x, mu_y) : (int, int)
        Center along X, Y
    (sigma_x, sigma_y) : (float, float)
        Standard deviation along X, Y

    Return
    ------
    A 2D mesh of topography with size (nx, ny)
    """
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    if mu_x == None:
        mu_x = nx // 2
    if mu_y == None:
        mu_y = ny // 2
    if sigma_x == None:
        sigma_x = nx / 8
    if sigma_y == None:
        sigma_y = ny / 8
    return np.exp(-((x - mu_x) ** 2 / (2 * sigma_x ** 2)
                  + (y - mu_y) ** 2 / (2 * sigma_x ** 2)))


# Make directories
Path('input').mkdir(parents=True, exist_ok=True)
Path('output_ckp').mkdir(parents=True, exist_ok=True)
Path('output_sfc').mkdir(parents=True, exist_ok=True)
Path('output_vol').mkdir(parents=True, exist_ok=True)

# Parameters of simulation
tmax = 5.0
dh = 45.0  # bottom dh
dt = 0.001
use_topography = 1

# Create mesh
# Top block
nx_0 = 360
ny_0 = 240
nz_0 = 64
vp = 4000 * np.ones((nz_0, ny_0, nx_0), dtype='float32')
vs = 2500 * np.ones((nz_0, ny_0, nx_0), dtype='float32')
rho = 2500 * np.ones((nz_0, ny_0, nx_0), dtype='float32')
vs[:60, :, :] = 250 
vp[:60, :, :] = 400
rho[:60, :, :] = 1800 
np.stack([vp, vs, rho], axis=-1).tofile('mesh_0')

# Bottom block
assert (nx_0 % 3 == 0) and (ny_0 % 3 == 0)
nx_1 = nx_0 // 3
ny_1 = ny_0 // 3
nz_1 = 32
vp = 4000 * np.ones((nz_1, ny_1, nx_1), dtype='float32')
vs = 2500 * np.ones((nz_1, ny_1, nx_1), dtype='float32')
rho = 2500 * np.ones((nz_1, ny_1, nx_1), dtype='float32')
np.stack([vp, vs, rho], axis=-1).tofile('mesh_1')


# Create source
srcx, srcy, srcz = nx_0 // 2, ny_0 // 3, nz_0 // 3
xoff, yoff, zoff = 0, 0, 0  # Better do not change
nsrc = 2
M0 = 1e14
nt_src = 5000
dt_src = dt  # Sometimes, the source dt is different from the simulation dt
dh_src = dh // 3  # Normally use top block spacing to define source location
t = np.arange(nt_src) * dt_src
tau = 0.4
stf = generate_gp(tau, t) * M0
gbuf = 1000

dir_src = "input/source.txt"
dir_src_xx = "input/source_xx"
dir_src_yy = "input/source_yy"
dir_src_zz = "input/source_zz"
dir_src_xz = "input/source_xz"
dir_src_yz = "input/source_yz"
dir_src_xy = "input/source_xy"

# file handles
fxx = open(dir_src_xx, 'w')
fyy = open(dir_src_yy, 'w')
fzz = open(dir_src_zz, 'w')
fxz = open(dir_src_xz, 'w')
fyz = open(dir_src_yz, 'w')
fxy = open(dir_src_xy, 'w')
np.random.seed(0)
with open(f"{dir_src}", "w") as f_coord:
    f_coord.write(
            f"2.0.0\n"
            f"file=input/source\n"
            f"degree=3\n"
            f"steps={nt_src}\n"
            f"stride=1\n"
            f"gpu_buffer_size={gbuf}\n"
            f"cpu_buffer_size={nt_src // gbuf}\n"
            f"length={nsrc}\n\n"
            f"coordinates\n")
    for i in range(nsrc):
        f_coord.write(f"0 {(srcx + 0 - xoff) * dh_src} {(srcy - yoff) * dh_src} {(zoff - srcz) * dh_src}\n")
        for j, f in enumerate([fxx, fyy, fzz, fxz, fyz, fxy]):
            stf.astype('float32').tofile(f)
            #np.multiply(stf, np.random.rand() - 0.5).astype('float32').tofile(f)

for fid in [fxx, fyy, fzz, fxz, fyz, fxy]:
    fid.close()

with open("source_0", "wb") as fid:
    for i in range(nsrc):
        idx = np.array((srcx  - xoff + 500, srcy - yoff, srcz), dtype='int32')
        idx.tofile(fid)
        #np.multiply(stf, np.random.rand() - 0.5).astype('float32').tofile(fid)
        np.tile(stf, (6, 1)).T.astype('float32').tofile(fid)

# Write topography file
if use_topography:
    pad = 8  # padding length, DO NOT CHANGE
    height = 50.0  # Maximum height of topography
    topo = height * gen_2d_gauss(nx_0 + pad * 2, ny_0 + pad * 2)
    header = np.array([nx_0, ny_0, pad], dtype='int32')
    with open('topography.bin', 'wb') as ftopo:
        header.tofile(ftopo)
        topo.T.astype('float32').tofile(ftopo)


# Create configuration
ND = 0  # ABC layers
ntiskp = 1
write_step = 100
newline = "\n"
with open('param.sh', 'w') as fid:
    fid.write(f' -X {nx_1} -Y {ny_1} -Z {nz_0},{nz_1} -x 3 -y 2 -G 2\n'
              f' --TMAX {tmax} --DH {dh} --DT {dt} --ND {ND} --ARBC 0.95\n'
              f' --IDYNA 0 --NSRC 2,0 --INSRC source --READ_STEP {nt_src}\n'
              f' --IFAULT 1 --NST {nt_src}\n'
              f' --MEDIASTART 2 --NTISKP {ntiskp} --WRITE_STEP {write_step}\n'
              f' --INVEL mesh --NVE 1 --NVAR 3\n'
              f' --NBGX 1,1 --NEDX {nx_0},{nx_1} --NSKPX 2,2\n'
              f' --NBGY 1,1 --NEDY {ny_0},{ny_1} --NSKPY 2,2\n'
              f' --NBGZ 1,1 --NEDZ {nz_0},{nz_1} --NSKPZ 2,2\n'
              f' -c output_ckp/ckp -o output_sfc\n'
              #f' --SOURCEFILE input/source.txt\n'
              #f' --RECVFILE input/receiver.txt\n'
              f'{" --INTOPO topography.bin{newline}" if use_topography else ""}'
              f' --FAC 1.0 --Q0 150. --EX 0.0 --FP 1.0')
    

# Create receivers 
# This part is used for within-motion project only
x0, y0 = nx_0 // 2, ny_0 // 2
width = 1
length = 1
depth = 2
dh = dh // 3
nsite = 2 * depth * (width + length - 2) + width * length if min(width, length) > 1 else width * length * (depth + 1)

# Hardcoded parameters
version = '2.0.0'
file_prefix = 'output_vol/receivers'
nsite = 10800  # Number of sites
steps = int(tmax / dt)
stride = ntiskp  # tiem steps to skip when writing
degree = 3  # interpolation
gpu_buffer_size = 50
cpu_buffer_size = 1  # Make sure steps % (gpu_buffer_size * cpu_buffer_size) == 0
num_writes = 1

assert steps % (gpu_buffer_size * cpu_buffer_size) == 0

# Write receiver.txt file
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
    count = 0
    for iz in range(depth):
        for iy in range(0, ny_0, 4):
            for ix in range(0, nx_0, 4):
                f_coord.write(f"0 {(ix - xoff) * dh} {(iy - yoff) * dh} {(zoff - iz) * dh}\n")
                count += 1 
    print(f'Count = {count}, modify the receiver.txt')
    """ 
    for iz in range(depth + 1):
        # output the whole volumn
        if iz == depth or min(length, width) == 1:
            for i in range(width):
                for j in range(length):
                    f_coord.write(f"0 {(x0 + j - xoff) * dh} {(y0 + i - yoff) * dh} {((zoff - iz) * dh - 10) * (1 if iz < 40 else 1)}\n")
            continue
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
    """

assert steps % (cpu_buffer_size * gpu_buffer_size) == 0
print("Done generating receivers.\n")


# Create submission script
with open('run.lsf', 'w') as fid:
    fid.write(f'#!/bin/bash\n'
              f'# Begin LSF Directives\n'
              f'#BSUB -J toy\n'
              f'#BSUB -P geo112\n'
              f'#BSUB -W 00:02\n'
              f'#BSUB -nnodes 1\n'
              f'#BSUB -alloc_flags "gpumps"\n'
              f'#BSUB -o run_%J.out\n'
              f'#BSUB -e run_%J.err\n'
              f'##BSUB -N\n\n'
  
              f'export OMP_NUM_THREADS=1\n'
              f'cd $LS_SUBCWD\n\n'
  
              f'module unload darshan-runtime\n'
              f'module load cuda\n'
              f'cat $0\n'
              f'args=$(cat param.sh)\n'
              f'echo $args\n\n'
  
              f'jsrun -n 6 -a 3 -c 3 -g 1 -r 6 -d cyclic /ccs/home/hzfmer/scratch/awp/pmcl3d $args')
              
              
