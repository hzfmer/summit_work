from awp_processing import awp

C = awp.Scenario("")
#fig=C.plot_slice_z(t = 10, iz=0, vmin=0.001, vmax=0.1, sx=srcidx[0], sy=srcidx[1], comp='X')
#sites = np.genfromtxt('stat.txt', skip_header=1)
v = C.read_syns_fromfile(fname='stat.txt', comps='xyz', chunk_size=20)
