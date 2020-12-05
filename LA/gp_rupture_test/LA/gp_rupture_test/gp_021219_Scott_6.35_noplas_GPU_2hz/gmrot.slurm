#!/bin/bash -l
#SBATCH -A geo112
#SBATCH -N 27
#SBATCH -t 0:05:00
#SBATCH -J gmrot


cd $SLURM_SUBMIT_DIR
date
echo "Running gmrot.pbs\n"
cat $0
srun -n 420 ~/file_back/programs/awp_tools/src/gmrot_rhea/gmrot_mpi
date
