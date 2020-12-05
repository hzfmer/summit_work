import numpy as np
import os
import sys
import pickle
import subprocess
import collections
from scipy.signal import resample


def resize(vel, vel_ref):
    length = len(vel_ref['t'])
    for comp in 'XYZ':
        vel[comp] = resample(vel[comp][vel['t'] <= vel_ref['t'][-1]], length)

def prepare_tf_misfit(model, vel_syn=None, fmin=0.15, fmax=5, exec_path='results1', IS_S2_REFERENCE='true', LOCAL_NORM='false'):
    if not vel_syn:
        with open('../results/vel_syn.pickle', 'rb') as fid:
            vel_syn = pickle.load(fid)

    try:
        misfit = dict()
        for j, k in enumerate(vel_syn[model].keys()):
            print(f"model={model}, site={k}")
            v = np.vstack((vel_syn['rec'][k]['t'], vel_syn[model][k]['X'], vel_syn[model][k]['Y'], vel_syn[model][k]['Z']))
            if np.max(v) < 1e-8:
                print("Error in velocities!")
                sys.exit(-1)
            np.savetxt(f'syn_{k}.dat', v.T, fmt='%.8f', delimiter=' ')
            with open('HF_TF-MISFIT_GOF', 'w') as fid:
                fid.write(f'{len(vel_syn["rec"][k]["X"])}\n{vel_syn["rec"][k]["dt"]}\n{fmin} {fmax}\n'
                          f'syn_{k}.dat\nrec_{k}.dat\n3\n.{IS_S2_REFERENCE}.\n.{LOCAL_NORM}.')
            command = ['./tf_misfits_gof', 'HF_TF-MISFIT_GOF']
            p = subprocess.run(command, capture_output=True, check=True)
            misfit[k] = np.genfromtxt('MISFIT-GOF.DAT')
            print(f"Done {j} / {len(vel_syn['rec'].keys())}")
            os.remove('MISFIT-GOF.DAT')
    except OSError as e:
        print("Error", e)
    finally:
        print("\n")
    return misfit


models = ['q100f06_orig_vs500', 'q100f00_orig_vs500', 'topo_q100f00_orig_vs500']
os.chdir('/ccs/home/hzfmer/scratch/high_f/la_habra_large_gpu_abc50/results1')

with open('../results/vel_syn.pickle', 'rb') as fid:
    vel_syn = pickle.load(fid)
try:
    with open('tf_misfit.pickle', 'rb') as fid:
        tf_misfit = pickle.load(fid)
except:
    tf_misfit = collections.defaultdict(dict)

freqs = [(0.15, 5), (0.15, 2.5), (2.5, 5)]

for i, model in enumerate(models):
    print(model)
    for k in vel_syn[model].keys():
        # k = sites 
        resize(vel_syn[model][k], vel_syn['rec'][k])
    for f in freqs:
        print(f)
        # if model in tf_misfit[f]:
        #     continue
        tf_misfit[f][model] = prepare_tf_misfit(model, vel_syn=vel_syn, fmin=f[0], fmax=f[1])
        # assert len(tf_misfit[f][model]) == len(vel_syn[model].keys())
        #print(tf_misfit.keys())

with open('tf_misfit.pickle', 'wb') as fid:
    pickle.dump(tf_misfit, fid, protocol=pickle.HIGHEST_PROTOCOL)

