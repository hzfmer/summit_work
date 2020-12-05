#!/usr/bin/env python3

import numpy as np
import pandas as pd

mx, my = 5940, 4860
nd = 250  # 2 km from the boundary

sites = pd.read_csv('la_habra_medium_statlist_5940.idx', sep=' ', header=None, index_col=False)


within = {}
j = 0
for i in range(len(sites)):
    if nd < sites.iloc[i, 1] < mx - nd and nd < sites.iloc[i, 2] < my - nd:
        #within[j] = [sites.iloc[i, 0], sites.iloc[i, 1], sites.iloc[i, 2], 1]
        within[j] = [sites.iloc[i, 1], sites.iloc[i, 2], 1]
        j += 1
seen = set()
print(type(within[0]), within[0])
with open('stat_medium.txt', 'w') as fid:
    fid.write(f'{len(within)}\n')
    for i in within:
        fid.write(" ".join(map(str, within[i])) + "\n")
        if (within[i][0], within[i][1]) in seen:
            print(f'Site {i} repeated', within[i][0], within[i][1])
        else:
            seen.add((within[i][0], within[i][1]))



