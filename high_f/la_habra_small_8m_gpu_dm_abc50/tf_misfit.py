import numpy as np
import os
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

models = ["noqf_orig", "qf06_orig", "noqf_orig_vs500", "noqf_s05h005l100"]
models = ["qf06_s05h005l100", "topo_noqf_orig_qs50", "noqf_orig_qs50", "topo_qf06_s05h005l100"]
models = ['topo_noqf_orig', 'topo_qf06_orig', "topo_qf06_s05h005l100_vs500"]

for dhyp, seed in zip([1,1], [1848640878, 387100462]):
    for case in ['noqf_orig', 'topo_noqf_orig', 'noqf_s05h005l100', 'topo_qf06_s05h005l100']:
        models = [f'dhyp{dhyp:.2f}_s{seed}_{case}'] + models

for dhyp, seed in zip([2,2], [372823598, 462574446]):
    for case in ['noqf_orig', 'topo_noqf_orig', 'noqf_s05h005l100', 'topo_qf06_s05h005l100']:
        models = [f'dhyp{dhyp:.2f}_s{seed}_{case}'] + models
# for k in vel_syn['rec'].keys():
#     v = np.vstack((np.arange(0, tpad, dt, dtype='float32'), vel_syn['rec'][k]['X'],
#                    vel_syn['rec'][k]['Y'], vel_syn['rec'][k]['Z']))
#     print(v.shape)
#     np.savetxt(f'results/rec_{k}.dat', v.T, fmt='%.8f', delimiter=' ')

os.chdir('/ccs/home/hzfmer/scratch/high_f/la_habra_small_8m_gpu_dm_abc50/results1')

with open('../results/vel_syn.pickle', 'rb') as fid:
    vel_syn = pickle.load(fid)
tf_misfit = collections.defaultdict(dict)
    
freqs = [(0.15, 5), (0.15, 2.5), (2.5, 5)]
tpad, dt = 35, 0.02

for i, model in enumerate(models):
    print(model)
    for k in vel_syn[model].keys():
        v = np.vstack((np.arange(0, tpad, dt, dtype='float32'), vel_syn[model][k]['X'],
                        vel_syn[model][k]['Y'], vel_syn[model][k]['Z']))
        np.savetxt(f'syn_{k}.dat', v.T, fmt='%.8f', delimiter=' ')           
    for f in freqs:
        print(f)
        tf_misfit[f][model] = prepare_tf_misfit(model, vel_syn=vel_syn, fmin=f[0], fmax=f[1])
        #print(tf_misfit.keys())

with open('tf_misfit.pickle', 'wb') as fid:
    pickle.dump(tf_misfit, fid, protocol=pickle.HIGHEST_PROTOCOL)
