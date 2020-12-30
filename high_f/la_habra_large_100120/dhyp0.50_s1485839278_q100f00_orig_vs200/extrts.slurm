#!/bin/bash
#SBATCH -J extrts       # Job name
#SBATCH -o extrts_%j.out          # Name of stdout output file
#SBATCH -e extrts_%j.err          # Name of stderr error file
#SBATCH -N 9              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 00:25:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

# Require IN3D.out
echo "Running extrts, Starting at `date`"
cat $0
srun -n 259  /ccs/home/hzfmer/scratch/awp_tools/bin_andes/extrts
echo "Ending at `date`"
