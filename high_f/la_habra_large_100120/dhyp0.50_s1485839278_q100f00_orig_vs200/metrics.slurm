#!/bin/bash
#SBATCH -J metrics       # Job name
#SBATCH -o metrics_%j.out          # Name of stdout output file
#SBATCH -e metrics_%j.err          # Name of stderr error file
#SBATCH -N 10              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 01:40:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo "Running metrics, Starting at `date`"
cat $0
srun -n 320  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 1 -n 6 -g
#srun -n 320  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 2.5 -n 6
# srun -n 320  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 2.5 -h 5 -n 6
#srun -n 320  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 5 -n 6

echo "Ending at `date`"
