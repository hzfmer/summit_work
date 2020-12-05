#!/usr/bin/env python

import os
import matplotlib
import numpy as np
import struct

file_mmt = "momrate_LaHabra.zf.100m.rev.bin"
file_out = "mesh_src.dat"
fidout = open(file_out,'w')
fidin = open(file_mmt, 'rb')

nsrc = 1560
nt = 5000
nmom = 6
nseek = 4 * nmom * nt

for i in range(nsrc):
    dat = struct.unpack('3I', fidin.read(12))
#   np.savetxt(file_out, dat)
    fidout.write(" ".join(str(x) for x in dat) + "\n") 
    fidin.seek(nseek,1)
       
fidin.close()
fidout.close()
    


