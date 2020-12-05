#!/usr/bin/env python
'''
Gaussian source description from:
https://books.google.com/books?id=F-3LDAAAQBAJ&pg=PA295&lpg=PA295&dq=seismic+gaussian+point+source&source=bl&ots=TYLyC_7xih&sig=HwLOv9h3erTyOSsYdjxdjpR5r2Y&hl=en&sa=X&ved=0ahUKEwjhnffx2-HbAhXgJDQIHRG7Aj84ChDoAQg_MAk#v=onepage&q&f=false

Center frequency: 4Hz
Time shift: 0.5s
'''

import numpy as np

dt = 0.0004
t = np.arange(0, 1+dt, dt)
T = 0.25
t0 = 0.5
p = 2 * np.pi * ((t - t0) / T)
A = 10**19
s = np.array([601, 238, 601], dtype='i')

F = np.array(-np.exp(-p ** 2 / np.pi ** 2) * np.sin(p), dtype='f')
F = np.tile(F, (3,1)).T
F = np.hstack((F, np.zeros_like(F))).flatten()
print(F.dtype)
src = open('source_0', 'wb')
src.write(s.tostring())
F.tofile(src)

src.close()
