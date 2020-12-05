#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfilt, sosfiltfilt
from scipy.fftpack import fft
import os

plt.tight_layout()
rcparams = {'font.size': 16,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 14,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'lines.linewidth': 2}
plt.rcParams.update(rcparams)


def filt_B(data_in, fs, lowcut=0, highcut=1, order=4, causal=True):
    if highcut >= fs/2:
        return data_in
    sz = data_in.shape
    if len(sz) > 1 and sz[0] > sz[1]:
        data_in = data_in.T
    data_out = np.zeros_like(data_in)
    nyq = fs / 2
    low = lowcut / nyq
    high = highcut / nyq
    if low == 0:  # lowpass
        sos = butter(order, high, analog=False, btype='low', output='sos')
    else:   # bandpass
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    if len(sz) > 1:
        for i in range(data_in.shape[0]):
            if causal == True:
                data_out[i, :] = sosfilt(sos, data_in[i, :])
                data_out[i, :] = sosfilt(sos, data_out[i, :])
            else:
                data_out[i, :] = sosfiltfilt(sos, data_in[i, :])
        if sz[0] > sz[1]:
            data_out = data_out.T
    else:
        if causal == True:
            data_out = sosfilt(sos, data_in)
            data_out = sosfilt(sos, data_out)
        else:
            data_out = sosfiltfilt(sos, data_in)
    return data_out

def rotate_v(vx, vy, theta):
    theta = np.deg2rad(theta)
    vx_out = vx * np.cos(theta) - vy * np.sin(theta)
    vy_out = vx * np.sin(theta) + vy * np.cos(theta)
    return vx_out, vy_out

def set_plot(v_cs, vx_awp, vy_awp, Mag, site, ix, iy, dist):
    vx_awp = filt_B(vx_awp, 20, 0.01, 1)
    vy_awp = filt_B(vy_awp, 20, 0.01, 1)
    vx_awp, vy_awp = rotate_v(vx_awp, vy_awp, -35)

    for i in range(2):
        vy_cs = filt_B(v_cs[0, :], 20, 0.01, 1)
        vx_cs = filt_B(v_cs[1, :], 20, 0.01, 1)

    fig, ax = plt.subplots(2, 1, figsize=(6,6))
    ax[0].plot(dt_cs * np.arange(len(vx_cs)), vx_cs, 'k', linewidth=1.5, label="CS")
    ax[0].plot(dt_awp * np.arange(len(vx_awp)), vx_awp, 'r--', dashes=(5, 1),  linewidth = 1, label="AWP")
    ax[0].set_xlabel('Time(s)')
    ax[0].set_ylabel('VX(m/s)')
    ax[0].set_title(f"M={Mag}, {site}, DIST={dist:.1f}km")
    ax[0].legend()

    ax[1].plot(dt_cs * np.arange(len(vy_cs)), vy_cs, 'k', linewidth=1.5, label="CS")
    ax[1].plot(dt_awp * np.arange(len(vy_awp)), vy_awp,  'r--', dashes=(5, 1), linewidth = 1, label="AWP")
    ax[1].set_xlabel('Time (s)')
    ax[1].set_ylabel('VY (m/s)')
    ax[1].legend()
    plt.subplots_adjust(hspace=0.3)

    fig.savefig(f'../results/{M}_{site}_{ix}_{iy}_rot_plas_2hz.png', dpi=400, bbox_inches='tight', 
            pad_inches=0.1)
    
    fvx_awp = np.abs(fft(vx_awp)) * dt_awp
    fvy_awp = np.abs(fft(vy_awp)) * dt_awp
    fvx_cs = np.abs(fft(vx_cs)) * dt_cs
    fvy_cs = np.abs(fft(vy_cs)) * dt_cs
    f_awp = np.linspace(0, 1 / 2 / dt_awp, len(vx_awp) // 2)
    f_cs = np.linspace(0, 1 / 2 / dt_cs, len(vx_cs) // 2)

    fig, ax = plt.subplots(2, 1, figsize=(6,6))
    ax[0].plot(f_cs, fvx_cs[0:len(fvx_cs)//2], 'k', linewidth=1.5, label="CS")
    ax[0].plot(f_awp, fvx_awp[0:len(fvx_awp)//2], 'r--', dashes=(5, 1), linewidth = 1, label="AWP")
    ax[0].set_xlabel('Frequency(Hz)')
    ax[0].set_ylabel('Amplitude')
    ax[0].set_title(f"M={Mag}, {site}, DIST={dist:.1f}km")
    ax[0].legend()
    ax[0].set_xlim((0, 1))

    ax[1].plot(f_cs, fvy_cs[0:len(fvy_cs)//2],'k', linewidth=1.5, label="CS")
    ax[1].plot(f_awp, fvy_awp[0:len(fvy_awp)//2], 'r--', dashes=(5, 1), linewidth = 1, label="AWP")
    ax[1].set_xlabel('Frequency (Hz)')
    ax[1].set_ylabel('Amplitude')
    ax[1].legend()

    plt.subplots_adjust(hspace=0.3)
    plt.xlim((0, 1))
    fig.savefig(f'../results/fft_{M}_{site}_{ix}_{iy}_rot_plas_2hz.png', dpi=400, bbox_inches='tight', 
            pad_inches=0.1)
    plt.close('all')
    return True

nt6, nx6, nz6 = 1000, 354, 212
nt7, nx7, nz7 = 2000, 980, 240
nt8, nx8, nz8 = 4000, 5116, 220
dt_cs = 0.05
dt_awp = 0.005 * 10
dh = 0.1

seis6 = np.load('../../data_for_Zhifeng/seis6.npy')[()]
seis7 = np.load('../../data_for_Zhifeng/seis7.npy')[()]
seis8 = np.load('../../data_for_Zhifeng/seis8.npy')[()]

f = open('../../data_for_Zhifeng/site_names6.txt', 'r')
stname6 = f.readlines()
stname6 = [i.strip('\n') for i in stname6] 
f.close()

f = open('../../data_for_Zhifeng/site_names7.txt', 'r')
stname7 = f.readlines()
stname7 = [i.strip('\n') for i in stname7] 
f.close()

f = open('../../data_for_Zhifeng/site_names8.txt', 'r')
stname8 = f.readlines()
stname8 = [i.strip('\n') for i in stname8] 
f.close()

f_lonlat6 = open('../gp_021219_Scott_6.35_plas_2hz/stat.txt','r')
f_lonlat7 = open('../gp_021219_Scott_7.35_plas_2hz/stat.txt','r')
f_lonlat8 = open('../gp_021219_Scott_8.45_plas_2hz/stat.txt','r')
nsite6 = f_lonlat6.readline()
ll6 = f_lonlat6.readlines()
ll6 = [i.strip('\n').split(' ')[:-1] for i in ll6]
ll6 = [[np.int(x) for x in l] for l in ll6]
nsite7 = f_lonlat7.readline()
ll7 = f_lonlat7.readlines()
ll7 = [i.strip('\n').split(' ')[:-1] for i in ll7]
ll7 = [[np.int(x) for x in l] for l in ll7]
nsite8 = f_lonlat8.readline()
ll8 = f_lonlat8.readlines()
ll8 = [i.strip('\n').split(' ')[:-1] for i in ll8]
ll8 = [[np.int(x) for x in l] for l in ll8]

dist6 = np.fromfile('site_rjb_M6.35.bin', dtype='float32')
dist7 = np.fromfile('site_rjb_M7.35.bin', dtype='float32')
dist8 = np.fromfile('site_rjb_M8.45.bin', dtype='float32')


sgms = os.listdir('../../data_for_Zhifeng/seismograms')
for i, sgm in enumerate(sgms):
    site = sgm.split('_')[1]

    if sgm.split('_')[2] == '68':
        M = 8.45
        idx = stname8.index(site)
        cs = seis8[site] / 100 
        ix = np.int(ll8[idx][0])
        iy = np.int(ll8[idx][1])
        s_dist = dist8[idx]
        try:
            awp_x = np.loadtxt(f'../gp_021219_Scott_8.45_plas_2hz/output_sfc/SX{ix:04d}_{iy:04d}_0001.dat')
            awp_y = np.loadtxt(f'../gp_021219_Scott_8.45_plas_2hz/output_sfc/SY{ix:04d}_{iy:04d}_0001.dat')
        except: 
            continue
        set_plot(cs[:nt8], awp_x, awp_y, '8.45', site, ix, iy, s_dist)
        print(f"\r{i}: Plotted {M}_{site}_{ix}_{iy}_2hz", end=' ', flush=True)

    if sgm.split('_')[2] == '76':
        M = 6.35
        idx = stname6.index(site)
        cs = seis6[site] / 100
        ix = np.int(ll6[idx][0])
        iy = np.int(ll6[idx][1])
        s_dist = dist6[idx]
        try:
            awp_x = np.loadtxt(f'../gp_021219_Scott_6.35_plas_2hz/output_sfc/SX{ix:04d}_{iy:04d}_0001.dat')
            awp_y = np.loadtxt(f'../gp_021219_Scott_6.35_plas_2hz/output_sfc/SY{ix:04d}_{iy:04d}_0001.dat')
        except:
            continue
        set_plot(cs[:nt6], awp_x, awp_y, '6.35', site, ix, iy, s_dist)
        print(f"\r{i}: Plotted {M}_{site}_{ix}_{iy}_2hz", end=' ', flush=True)
    
    if sgm.split('_')[2] == '128':
        M = 7.35
        idx = stname7.index(site)
        cs = seis7[site] / 100 
        ix = np.int(ll7[idx][0])
        iy = np.int(ll7[idx][1])
        s_dist = dist7[idx]
        try:
            awp_x = np.loadtxt(f'../gp_021219_Scott_7.35_plas_2hz/output_sfc/SX{ix:04d}_{iy:04d}_0001.dat')
            awp_y = np.loadtxt(f'../gp_021219_Scott_7.35_plas_2hz/output_sfc/SY{ix:04d}_{iy:04d}_0001.dat')
        except:
            continue
        set_plot(cs[:nt7], awp_x, awp_y, '7.35', site, ix, iy, s_dist)
        print(f"\r{i}: Plotted {M}_{site}_{ix}_{iy}_2hz", end=' ', flush=True)
    plt.close('all')

