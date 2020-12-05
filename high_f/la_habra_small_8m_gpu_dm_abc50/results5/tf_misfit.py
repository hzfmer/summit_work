import numpy as np
import os
import sys
import pickle
import subprocess
import collections


def prepare_tf_misfit(model, vel_syn=None, fmin=0.15, fmax=5, exec_path='results1', IS_S2_REFERENCE='true', LOCAL_NORM='false'):
    if not vel_syn:
        with open('../results/vel_syn.pickle', 'rb') as fid:
            vel_syn = pickle.load(fid)

    try:
        tf_misfit = dict()
        for j, k in enumerate(vel_syn['rec'].keys()):
            with open('HF_TF-MISFIT_GOF', 'w') as fid:
                fid.write(f'{len(vel_syn[model][k]["X"])}\n{vel_syn[model][k]["dt"]}\n{fmin} {fmax}\n'
                          f'syn_{k}.dat\nrec_{k}.dat\n3\n.{IS_S2_REFERENCE}.\n.{LOCAL_NORM}.')
            command = ['./tf_misfits_gof', 'HF_TF-MISFIT_GOF']
            p = subprocess.call(command,
                                stdout = subprocess.PIPE,
                                stdin  = subprocess.PIPE,
                                stderr = subprocess.STDOUT )
            tf_misfit[k] = np.genfromtxt('MISFIT-GOF.DAT')
            print(f"\rDone {j} / {len(vel_syn['rec'].keys())}", end='\r', flush=True)
            os.remove('MISFIT-GOF.DAT')
    except OSError as e:
        print("Error", e)
    finally:
        print("\n")
    return tf_misfit

models = [c + "_s05h005l100" for c in ["topo_q50f02", "topo_q50f04", "topo_q50f06", "topo_q50f08", "topo_q100f00", "topo_q100f02", "topo_q100f04", "topo_q100f08"]]
# for k in vel_syn['rec'].keys():
#     v = np.vstack((np.arange(0, tpad, dt, dtype='float32'), vel_syn['rec'][k]['X'],
#                    vel_syn['rec'][k]['Y'], vel_syn['rec'][k]['Z']))
#     print(v.shape)
#     np.savetxt(f'results/rec_{k}.dat', v.T, fmt='%.8f', delimiter=' ')

os.chdir('/ccs/home/hzfmer/scratch/high_f/la_habra_small_8m_gpu_dm_abc50/results5')

with open('../results/vel_syn.pickle', 'rb') as fid:
    vel_syn = pickle.load(fid)
try:
    with open('tf_misfit.pickle', 'rb') as fid:
        tf_misfit = pickle.load(fid)
except:
    tf_misfit = collections.defaultdict(dict)

freqs = [(0.15, 5), (0.15, 2.5), (2.5, 5)]
tpad, dt = 35, 0.02

for i, model in enumerate(models):
    print(model)
    for k in vel_syn[model].keys():
        v = np.vstack((np.arange(0, tpad, dt, dtype='float32'), vel_syn[model][k]['X'],
                        vel_syn[model][k]['Y'], vel_syn[model][k]['Z']))
        if np.max(v) < 1e-8:
            print("Error in velocities!")
            sys.exit(-1)
        np.savetxt(f'syn_{k}.dat', v.T, fmt='%.8f', delimiter=' ')
    for f in freqs:
        print(f)
        #if model not in tf_misfit[f]:
        tf_misfit[f][model] = prepare_tf_misfit(model, vel_syn=vel_syn, fmin=f[0], fmax=f[1])
        #print(tf_misfit.keys())

with open('tf_misfit.pickle', 'wb') as fid:
    pickle.dump(tf_misfit, fid, protocol=pickle.HIGHEST_PROTOCOL)

