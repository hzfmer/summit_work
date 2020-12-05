#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys

nt = 100
dt = 0.005
fname = sys.argv[1]
rank = fname.split('_')[0]

f_tpsrc = open(f'../tpsrc/{rank.replace("fault", "tpsrc")}', 'rb')
npsrc = np.frombuffer(f_tpsrc.read(4), dtype='i')[0]

print(npsrc)

data = np.fromfile(fname, dtype='f').reshape(6, npsrc, nt)
print(np.max(data))
idx = np.argmax(np.abs(data))
i = idx // (npsrc * nt)
j = (idx % (npsrc * nt)) // nt
fig, ax = plt.subplots(figsize=(6,6))
ax.plot(np.linspace(0, nt *  dt, nt), data[i, j, :])
print(np.max(data[i, j, :]))
fig.savefig(f"plot_max_subfaults.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

