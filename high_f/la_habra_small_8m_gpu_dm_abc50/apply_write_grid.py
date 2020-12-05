import numpy as np
import sys
import subprocess

nx, ny, nz, mz, NZ = 3500, 3500, 108, 180, 1000
NZ = (NZ - nz) // 4 // 3 * 4 * 3 + nz - 7  # 7 layers in burface block
mx, my = nx // 54 * 54, ny // 54 * 54
px, py = (nx - mx) // 2, (ny - my) // 2
pad = 8
length = 4 * 3 * nx * ny

try:
    # raise Exception
    fmesh = open('mesh.bin', 'rb')
    fout = open('mesh_1.bin','wb')
    fmesh.seek(3 * 4 * mx * my * (nz - 7))
    mesh = np.frombuffer(fmesh.read(3 * 4 * mx * my * 7), dtype='float32').reshape(7, my, mx, 3)
    mesh[::3, ::3, ::3].tofile(fout)
    subprocess.call('cat bottom_cvm.bin >> mesh_1.bin', shell=True)

except IOError:
    with open('mesh', 'rb') as fmesh,  \
        open('cvm.bin', 'w') as top_cvm,  \
        open('bottom_cvm.bin', 'w') as bot_cvm:
        for i in range(nz):
            buf = fmesh.read(length)
            cvm = np.frombuffer(buf, dtype='float32').reshape(ny, nx, 3)[py : py + my, px : px + mx, :]
            cvm.tofile(top_cvm)
        for i in range(nz, mz):
            buf = fmesh.read(length)
            cvm = np.frombuffer(buf, dtype='float32').reshape(ny, nx, 3)[py : py + my, px : px + mx, :]
            cvm.tofile(top_cvm)
            if i % 3 == 2:
                cvm[::3, ::3, :].tofile(bot_cvm)
        for i in range(mz, NZ):
            buf = fmesh.read(length)
            if i % 3 == 2:
                cvm = np.frombuffer(buf, dtype='float32').reshape(ny, nx, 3)[py : py + my, px : px + mx, :]
                cvm[::3, ::3, :].tofile(bot_cvm)
