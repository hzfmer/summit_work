import numpy as np
from Pathlib2 import Path 

fine_fname = "../cvm/la_habra_large_cvmsi_8m.media"
fine_add_fname = "../cvm/la_habra_large_cvmsi_8m_25.media"
coarse_fname = "../cvm/la_habra_large_cvmsi_8m_bot.media"

# Output filenames
top_fname = "mesh_large_8m_orig.bin_0"
center_fname = "mesh_large_8m_orig.bin_1"
bot_fname = "mesh_large_8m_orig.bin_2"

if ssh:
    # SSH: s05h005l100
    std = 10
    hurst = 5
    l = 100
    fine_fname = f"mesh_large_8m_s{std:02d}h{hurst:03d}l{l}.bin"
    top_fname = f"mesh_large_8m_s{std:02d}h{hurst:03d}l{l}.bin_0"
    center_fname = f"mesh_large_8m_s{std:02d}h{hurst:03d}l{l}.bin_1"
    write_bot = 0

write_top, write_center, write_bot = 1, 1, 0
nvar = 3
bitfloat = 4
NZ = 1250 
#=================================
# Generate the first block of mesh
# 0 - 1472m : from the finest 8m concatenated UCVM mesh (hom/het)
# Index : 0 - 183 
nx, ny, nz = 19440, 14904, 184
nlayer = nx * ny * nvar 
if write_top:
    with open(top_fname, 'wb') as fout:
        for i in range(nz):
            data = np.fromfile(fine_fname, dtype='float32', count=nlayer, 
                    offset=nlayer*bitfloat*i).reshape(ny, nx, nvar)
            data.tofile(fout)


#=================================
# Generate the second block of mesh
# 1216m - 9976m : from the finest 8m concatenated UCVM mesh (hom/het)
# Index: 152 - 1249
# 10000 - 10312m : from the finest 8m UCVM mesh **_8m_25.media
# Index: 0 - 39 
mx, my = nx // 3, ny // 3
mz = 380
if write_center:
    with open(center_fname, 'wb') as fout:
        for i in range(nz - 8, NZ, 3):
            print(f"Processing layer {i} / {NZ} in the center block", end="\r", flush=True)
            data = np.fromfile(fine_fname, dtype='float32', count=nlayer,
                    offset=nlayer*bitfloat*i).reshape(ny, nx, nvar)
            data[1::3, ::3, :].tofile(fout)
        for i in range(0, 40, 3):
            data = np.fromfile(fine_add_fname, dtype='float32', count=nlayer,
                    offset=nlayer*bitfloat*i).reshape(ny, nx, nvar)
            data[1::3, ::3, :].tofile(fout)


#=================================
# Generate the bot block of mesh 
# 10144m - 10288m : from the finest 8m UCVM mesh **_8m_25.media 
# Index: 18 - 36
# 10360 - 59032m : from the coarsest 72m UCVM mesh 
kx, ky = nx // 9, ny // 9 
klayer = kx * ky * nvar 
kz = 680  # the number of layers to append from the UCVM coarsest mesh 
if write_bot: 
    with open(bot_fname, 'wb') as fout: 
        for i in range(18, 36+1, 9): 
            data = np.fromfile(fine_add_fname, dtype='float32', count=nlayer,  offset=nlayer*bitfloat*i).reshape(ny, nx, nvar)
            data[4::9, ::9, :].tofile(fout)
        for i in range(kz - 3):  # 3 layers added above
            data = np.fromfile(coarse_fname, dtype='float32', count=klayer,
                    offset=klayer*bitfloat*i)
            data.tofile(fout) 

# Check continuity
if write_top or write_center: 
    if Path(top_fname).exists() and Path(center_fname).exists():
        check_mesh_cont(top_fname, center_fname, nx, ny, nz)

if write_center or write_bot:
    if Path(center_fname).exists() and Path(bot_fname).exists():
        check_mesh_cont(center_fname, bot_fname, mx, my, mz) 

