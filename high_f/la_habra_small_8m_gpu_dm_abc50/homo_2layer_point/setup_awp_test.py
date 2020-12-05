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

import numpy as np
from pathlib import Path
import sys
from scipy.signal import resample


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
tmax = 8.0
dh = 45.0  # bottom dh
dt = 0.001

if len(sys.argv) == 1:
    NPX = 2
    NPY = 2
else:
    NPX = int(sys.argv[1])
    NPY = int(sys.argv[2])
    if max(NPX, NPY) > 3:
        print(f'For one node test, using <= 3 GPUs along each derection.\nAborting!')
        sys.exit(-1)

# Create mesh
# Top block
nx_0 = 720
ny_0 = 450
nz_0 = 128
vp = 4000 * np.ones((nz_0, ny_0, nx_0), dtype='float32')
vs = 2500 * np.ones((nz_0, ny_0, nx_0), dtype='float32')
rho = 2500 * np.ones((nz_0, ny_0, nx_0), dtype='float32')
np.stack([vp, vs, rho], axis=-1).tofile('mesh_0')

# Bottom block
assert (nx_0 % 3 == 0) and (ny_0 % 3 == 0)
nx_1 = nx_0 // 3
ny_1 = ny_0 // 3
nz_1 = 64
vp = 4000 * np.ones((nz_1, ny_1, nx_1), dtype='float32')
vs = 2500 * np.ones((nz_1, ny_1, nx_1), dtype='float32')
rho = 2500 * np.ones((nz_1, ny_1, nx_1), dtype='float32')
np.stack([vp, vs, rho], axis=-1).tofile('mesh_1')
assert (nx_1 % NPX == 0) and (ny_0 % NPY == 0)

# Create source
srcx, srcy, srcz = nx_0 // 2, ny_0 // 3, nz_0 // 3  # artificially selected
xoff, yoff, zoff = 0, 0, 0  # Better do not change
nsrc = 1
nt_src = 1000
M0 = 10 ** 15
dt_src = dt  # Sometimes, the source dt is different from the simulation dt
dh_src = dh // 3  # top block grid spacing
t = np.arange(nt_src) * dt_src
tau = 0.4
stf = M0 * generate_gp(tau, t)
gbuf = nt_src

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
    stf = resample_src(stf, int(nt_src * dt_src / dt), int(tmax / dt))
    print(stf.shape)
    for i in range(nsrc):
        f_coord.write(f"0 {(srcx + i - xoff) * dh_src} {(srcy - yoff) * dh_src} {(zoff - srcz) * dh_src}\n")
        for j, f in enumerate([fxx, fyy, fzz, fxz, fyz, fxy]):
            np.multiply(stf, np.random.rand() - 0.5).astype('float32').tofile(f)

for fid in [fxx, fyy, fzz, fxz, fyz, fxy]:
    fid.close()


# Write topography file
pad = 8  # padding length, DO NOT CHANGE
height = 50.0  # Maximum height of topography
topo = height * gen_2d_gauss(nx_0 + pad * 2, ny_0 + pad * 2)
topo = 30 * np.ones((nx_0 + pad * 2, ny_0 + pad * 2))
header = np.array([nx_0, ny_0, pad], dtype='int32')
with open('topography.bin', 'wb') as ftopo:
    header.tofile(ftopo)
    topo.T.astype('float32').tofile(ftopo)


# Create configuration
ND = 10  # ABC layers
ntiskp = 20
write_step = int(tmax / dt) // ntiskp // 2
for (NPX, NPY) in [(1,1), (2, 1), (1,3), (2, 3), (2,2), (3, 2), (3, 1)]:
    Path(f'output_sfc_{NPX}_{NPY}').mkdir(parents=True, exist_ok=True)
    with open(f'param_{NPX}_{NPY}.sh', 'w') as fid:
        fid.write(f' -X {nx_1} -Y {ny_1} -Z {nz_0},{nz_1} -x {NPX} -y {NPY} -G 2\n'
                f' --TMAX {tmax} --DH {dh} --DT {dt} --ND {ND} --ARBC 0.95\n'
                f' --IDYNA 0 --NSRC 0,0\n'
                f' --MEDIASTART 2 --NTISKP {ntiskp} --WRITE_STEP {write_step}\n'
                f' --INVEL mesh --NVE 1 --NVAR 3\n'
                f' --NBGX 1,1 --NEDX {nx_0},{nx_1} --NSKPX 4,2\n'
                f' --NBGY 1,1 --NEDY {ny_0},{ny_1} --NSKPY 2,2\n'
                f' --NBGZ 1,1 --NEDZ {nz_0},{nz_1} --NSKPZ 2,2\n'
                f' -c output_ckp/ckp -o output_sfc_{NPX}_{NPY}\n'
                f' --SOURCEFILE input/source.txt\n'
                f' --INTOPO topography.bin\n'
                f' --FAC 1.0 --Q0 150. --EX 0.0 --FP 1.0')
        


    # Create submission script
    with open(f'run_{NPX}_{NPY}.lsf', 'w') as fid:
        fid.write(f'#!/bin/bash\n'
                f'# Begin LSF Directives\n'
                f'#BSUB -J toy\n'
                f'#BSUB -P geo112\n'
                f'#BSUB -W 00:{max(4, int(6 / NPX) * int(4 / NPY)):02d}\n'
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
                f'args=$(cat param_{NPX}_{NPY}.sh)\n'
                f'echo $args\n\n'
    
                f'jsrun -n {NPX * NPY} -a 3 -c 3 -g 1 -r {NPX *NPY} -d cyclic /ccs/home/hzfmer/awp_highf/pmcl3d $args')
                
                
