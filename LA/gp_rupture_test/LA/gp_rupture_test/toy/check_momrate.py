import numpy as np
import matplotlib.pyplot as plt

nx, nz, nt = 354, 212, 4000
dt = 0.005
ntskp = 10

slipr = np.fromfile('sliprate.bin', dtype='f').reshape(-1, nz, nx)
nt_slip = slipr.shape[0]
with open("momrate.dat", 'rb') as f_mom:
    fig, ax = plt.subplots(4, 1, figsize=(6, 6))
    for i in range(4):
        np.random.seed()
        ix = np.random.randint(nx)
        iz = np.random.randint(nz)
        disp = iz * nx + ix
        print(disp)
        f_mom.seek((12 + 4 * 6 * nt) * disp + 12, 0)
        # momrate has 6 components
        buf_mom = np.frombuffer(f_mom.read(4 * 6 * nt), dtype='f')[::6]
        max_mom = max(abs(buf_mom))
        sum_mom = buf_mom.sum() * dt
        buf_slip = slipr[:, iz, ix]
        max_slip = max(abs(buf_slip))
        if max_slip < 1e-20 or max_mom < 1e-20:
            continue
        buf_mom = buf_mom / max_mom
        buf_slip = buf_slip / max_slip
        print("moment, max_momrate, max_sliprate: ", sum_mom, max_mom, max_slip)
        ax[i].plot(np.arange(0, nt * dt, dt), buf_mom, 'k', linewidth=1.5)
        ax[i].plot(np.arange(0, nt_slip * dt * ntskp, dt * ntskp), buf_slip, 'r--', linewidth=2)
        ax[i].set_ylabel('Amplitude')
    ax[0].legend(['momrate', 'sliprate'], loc=1)
    ax[3].set_xlabel('Time (s)')
    fig.savefig(f"check_momrate.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
