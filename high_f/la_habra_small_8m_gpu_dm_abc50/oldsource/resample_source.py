#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
import numpy as np
import re
import fileinput

with open("param.sh") as fid:
    for line in fid:
        params = fid.readline()
        if "TMAX" in params:
            tmax = float(re.findall(r'(?<=TMAX )\d+\.?\d+(?= )', params)[0])
        if "DT" in params:
            dt_upsrc = float(re.findall(r'(?<=DT )\d+\.?\d+(?= )', params)[0])
            break
dt = 0.001
nt = 5000
tsrc = dt * nt
nsrc = 1560
npend = int((tmax - tsrc) / dt_upsrc) 
print(f"TMAX = {tmax}, tsrc={tsrc}, dt_upsrc={dt_upsrc}")
ntupsrc = int(tsrc / dt_upsrc)
zeros = np.zeros((nsrc, npend), dtype='float32')
upsrc = np.zeros((nsrc, ntupsrc), dtype='float32')
print(f"Upsample to {ntupsrc}, in total {npend + ntupsrc} time steps")

#with fileinput('input/source.txt', inplace=True) as fid:
gbuf = 1000
with open('../input/source.txt', 'r') as fid, open('input/source.txt', 'w') as fout:
    for line in fid:
        if "steps" in line:
            line = f"steps={npend + ntupsrc}\n"
        if "gpu" in line:
            gbuf = int(line.split('=')[1])
        if "cpu" in line:
            line = f"cpu_buffer_size={(npend + ntupsrc) // gbuf}\n"
        fout.write(line)

for comp in ["xx", "yy", "zz", "xz", "yz", "xy"]:
    try:
        src = np.fromfile(f'../input/source_{comp}', dtype='float32').reshape(nsrc, nt)
    except Exception as e:
        print("Error: \n", e)
    for i in range(nsrc):
        upsrc[i, :] = np.interp(np.arange(0, tsrc, dt_upsrc), np.arange(0, tsrc, dt), src[i, :])
    src = np.concatenate((upsrc, zeros), axis=1)
    if np.shape(src)[1] != ntupsrc + npend:
        print("Error! Mismatch in length of upsampled source time series")
        break
    src.tofile(f'input/source_{comp}')

