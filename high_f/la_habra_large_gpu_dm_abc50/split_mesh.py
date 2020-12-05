#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import sys

print("Usage: ./split_mesh.py nx ny nz NZ f_mesh f_p0 f_p1")
if len(sys.argv) < 8:
    print("Incorrect input\n")
    sys.exit(-1)
nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])
NZ = int(sys.argv[4])
f_mesh = sys.argv[5]
f_p0 = sys.argv[6]
f_p1 = sys.argv[7]
remainder = (nz - 8) % 3

# nx, ny, nz, NZ = 3456, 3456, 108, 1500
#NZ = (NZ - nz) // 4 // 3 * 4 * 3 + nz - 8  # 8 layers in surface block
length = 4 * 3 * nx * ny
print(f"Shape of the top block: ({nx}, {ny}, {nz})\n"
      f"The number of layers of the second block: {(NZ + 8 - nz) // 3}")

with open(f_mesh, 'rb') as fmesh, \
        open(f_p0, 'wb') as f0, \
        open(f_p1, 'wb') as f1:
            for i in range(nz - 8):
                buf = fmesh.read(length)
                cvm = np.frombuffer(buf, dtype='float32').reshape(ny, nx, 3)
                cvm.tofile(f0)
            for i in range(nz - 8, nz):
                buf = fmesh.read(length)
                cvm = np.frombuffer(buf, dtype='float32').reshape(ny, nx, 3)
                cvm.tofile(f0)
                if i % 3 == remainder:
                    cvm[1::3, ::3, :].tofile(f1)
            for i in range(nz, NZ):
                buf = fmesh.read(length)
                if i % 3 == remainder:
                    cvm = np.frombuffer(buf, dtype='float32').reshape(ny, nx, 3)
                    cvm[1::3, ::3, :].tofile(f1)


