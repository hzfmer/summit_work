import numpy as np

nbit = 4  # float bit
"""
nx, ny, nz, nvar = 19440, 14904, 160, 3
mx, my = nx * 2 // 3 + 1197, ny * 2 // 3
endx, endy = mx + 2592, my + 2160
nlayer = nx * ny * nvar
with open("mesh_0", "wb") as fid:
    for i in range(nz):
        data = np.fromfile("../mesh_large_8m_orig.bin_0", dtype='float32',
            count=nlayer, offset=nlayer * i * nbit).reshape(ny, nx, nvar)
        data[my : endy, mx : endx, :].tofile(fid)

nx, ny, nz, nvar = 19440 // 3, 14904 // 3, 380, 3
mx, my = nx * 2 // 3 + 1197 // 3, ny * 2 // 3
endx, endy = mx + 2592 // 3, my + 2160 // 3
nlayer = nx * ny * nvar
with open("mesh_1", "wb") as fid:
    for i in range(nz):
        data = np.fromfile("../mesh_large_8m_orig.bin_1", dtype='float32',
            count=nlayer, offset=nlayer * i * nbit).reshape(ny, nx, nvar)
        data[my : endy, mx : endx, :].tofile(fid)

nx, ny, nz, nvar = 19440 // 9, 14904 // 9, 60, 3
mx, my = nx * 2 // 3 + 1197 // 9, ny * 2 // 3
endx, endy = mx + 2592 // 9, my + 2160 // 9
nlayer = nx * ny * nvar
with open("mesh_2", "wb") as fid:
    for i in range(nz):
        data = np.fromfile("../mesh_large_8m_orig.bin_2", dtype='float32',
            count=nlayer, offset=nlayer * i * nbit).reshape(ny, nx, nvar)
        data[my : endy, mx : endx, :].tofile(fid)
"""
from awp_processing.check import check_mesh_cont
check_mesh_cont("mesh_0", "mesh_1", 2592, 2160, 160)
check_mesh_cont("mesh_1", "mesh_2", 2592 // 3, 2160 // 3, 380)
