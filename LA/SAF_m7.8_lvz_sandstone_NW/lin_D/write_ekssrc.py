#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from pylab import *
import os
from os.path import dirname as up
import struct
from scipy import signal
from sys import argv

def stf_transfer(stf, ratio, pw):
# stf: source time function
# ratio: control amplitude
# pw: power of exponential function, controlling the width.
#     6 is a preferrable value now, which makes it close to normal
#     distribution
    max_id = argmax(absolute(stf))
    stf_t = zeros_like(stf)
    if max_id < len(stf) - 1 and stf[max_id + 1] == 0:
        max_id = 0
    if max_id == 0:
        last_zr = 0
    else:
        try:
            # look for nearest zeros to the max_id
            # or just last_zr = argwhere(absolute(stf) > 0.01)[0]
            #         last_zr = argwhere(last_zr > max_id - 20), 20 can be repla            #                   ce by 1/dt
            last_zr = argwhere(
                            (isclose(stf[0:max_id],0,atol=max(5e-3,max(absolute(stf)) * 1e-2)))
                            &
                            (isclose(diff(stf[0:max_id+1]),0,atol=1e-3)))[-1]    
        except:
            last_zr = 0

    # ESSENTIAL float here
    # float is necessary, otherwise stupid old-fashioned py2.7 will 
    # get integers for compuations using rtime.
    rtime = float(max_id - last_zr)
    if rtime < 1:
        stf_t = stf.copy()
    else:
        t = arange(last_zr, len(stf))
        t_shift = (t - last_zr) / rtime
        factor = np.power(t_shift, pw)  * exp(- np.power(t_shift, pw))
        factor = 1 - factor * exp(1) * (1 - ratio)
        stf_t[t] = stf[t] * factor
    return stf_t

if len(argv) != 2:
    print("Usage: {0} rock".format(argv[0]))
    rock = "sandstone"
else:
    rock = argv[1]
nx = 2960
nz = 180
nt = 2000
dx = 0.1
dt = 0.06
pw = 6    # some magic value used in stf_transfer
          # 6 gives a narrower width of transfer function than 1
          # can be revised in future, perhaps distance-related
r_u = 0.95
r_w = 0.8
r0_u = 0.4
r0_w = 0.5
GSI = 50
if "sandstone" in rock:
    r_u = r_u
    r_w = r_w
    r0_u = 0.4
    r0_w = 0.5
    GSI = 50
elif "shale" in rock:
    r0_u = 0.3
    r0_w = 0.4
    GSI = 30

case = os.getcwd()[-1]
fid_stress = open(os.path.join(up(up(os.getcwd())),
                        'source_7.8_so_cvm_lvz/source_%s.txt' % case))
stress = float(fid_stress.readlines()[-1].strip('\n').split(':')[-1])

delta_r = (stress - 7e6) / 7e6 / 10
print(delta_r)
r_u = r_u - delta_r
r_w = r_w - delta_r

u_loc = GSI * 3 / 50; 
w_loc = u_loc + 0.5;

p_ratu = fromfile("output_dyn/maxrate_x_01.0Hz.bin", 'f').reshape(nz, nx)
p_ratw = fromfile("output_dyn/maxrate_z_01.0Hz.bin", 'f').reshape(nz, nx)
ratu = fromfile("output_dyn/ratu", 'f').reshape(nt, nz, nx)
ratw = fromfile("output_dyn/ratw", 'f').reshape(nt, nz, nx)

ratu_dep = flipud(mean(p_ratu, axis=1))
ratw_dep = flipud(mean(p_ratw, axis=1))

# looking for at least 5 smaller diff values
# pick 1/2 location as limit velocity 
#id_u = ceil(where(cumsum(absolute(diff(
#                ratu_dep)) < 0.01) == 2)[0][-1] / 2).astype(int)
#id_w = ceil(where(cumsum(absolute(diff(
#                ratw_dep)) < 0.01) == 2)[0][-1] / 2).astype(int)
id_w = argwhere(absolute(diff(ratw_dep)) < 0.03)
id_w = id_w[3]
id_u = id_w
lim_u = ratu_dep[round(asscalar(id_u/2))]
lim_w = ratw_dep[round(asscalar(id_w/2))]

# some magic random functions, who knows if work or not
# starting at about 0.45 on the surface, increasing to about 1(u)/0.9(w) (may
# affected by rake direction), adding some noises
# note that the depth column is stored with surface as the last.

# ratio on the surface, 0.45 for sandstone, 0.3 for shale
# to be modified according to rock strength 
#ratio_s = 0.45
#ratio_s = 0.3
ratio_u = r0_u + np.random.normal(0, 0.02, nz) + (r_u - r0_u) * (1 -
            exp(-linspace(0, u_loc * nz / id_u, nz)))
ratio_u = ratio_u[::-1]
ratio_w = r0_w + np.random.normal(0, 0.02, nz) + (r_w - r0_w) * (1 -
            exp(-linspace(0, w_loc * nz / id_w, nz)))
ratio_w = ratio_w[::-1]

fs = 1 / dt
fc = 1
[b, a] = signal.butter(4, fc / (fs / 2), 'low') 


ratu_t = zeros_like(ratu)
ratw_t = zeros_like(ratw)
for kk in range(nz):
    for ll in range(nx):
        tmp = signal.lfilter(b, a, ratu[:, kk, ll])
# can add a random function for each ll?
        if max(tmp) * ratio_u[kk] > lim_u and (nz-kk) <= id_u:
            ratu_t[:, kk, ll] = stf_transfer(tmp, lim_u / max(tmp), pw)
        else:
            ratu_t[:, kk, ll] = stf_transfer(tmp, ratio_u[kk], pw)

        tmp = signal.lfilter(b, a, ratw[:, kk, ll])
        if max(tmp) * ratio_w[kk] > lim_w and (nz-kk) <= id_w:
            ratw_t[:, kk, ll] = stf_transfer(tmp, lim_w / max(tmp), pw)
        else:
            ratw_t[:, kk, ll] = stf_transfer(tmp, ratio_w[kk], pw)

# output files, second replace for different rocks
output_u = open(os.path.join(os.getcwd(), "output_dyn/ratu").replace('lin','eks')
              .replace('sandstone', rock)
              ,'wb')
output_w = open(os.path.join(os.getcwd(), "output_dyn/ratw").replace('lin','eks')
              .replace('sandstone', rock)
              , 'wb')

for ii in range(nt):
    for jj in range(nz):
        output_u.write(struct.pack('<%df' % nx, *ratu_t[ii,jj,:]))
        output_w.write(struct.pack('<%df' % nx, *ratw_t[ii,jj,:]))
#output_u.write(struct.pack('<%df' % len(ratu_t), *ratu_t))
#output_w.write(struct.pack('<%df' % len(ratw_t), *ratw_t))

output_u.close()
output_w.close()


