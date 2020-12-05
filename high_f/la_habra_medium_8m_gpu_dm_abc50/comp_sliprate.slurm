#!/bin/bash
#SBATCH -J comp_sliprate       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 1              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 00:30:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

module load python/3.7.0-anaconda3-2018.12
echo "Running comp_sliprate, Starting at `date`"
srun -n 1 python comp_sliprate.py
echo "Ending at `date`"
