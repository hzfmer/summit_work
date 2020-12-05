#!/bin/bash
#SBATCH -J add_het_p0  # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 16              # Total # of nodes
##SBATCH -n 108           # Total # of mpi tasks
#SBATCH -t 00:05:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo "Running add_het_p0, Starting at `date`"
# srun -n 108 ./add_het_p0
# -n is number of layers
srun -n 500 ./add_het_p0
#srun -n 3072 --npernode 8 ./add_het_p0
echo "Ending at `date`"
