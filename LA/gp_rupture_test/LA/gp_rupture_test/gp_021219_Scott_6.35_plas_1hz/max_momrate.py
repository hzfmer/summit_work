#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from struct import pack
import os

npt = 2960 * 180
maxx = 1e8
with open("momrate.dat", 'rb') as f:
    for i in range(20):
        f.seek(12,1)
        for j in range(20):
            buffer = f.read(24000)
            tmp = np.frombuffer(buffer, dtype='f', count=6000)
            maxx = max([tmp.max(), maxx])
        print(np.abs(maxx))
    
