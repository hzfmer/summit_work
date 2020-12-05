#!/bin/bash
#SBATCH -J latlon2ij       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 1              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 00:20:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo "Running latlon2ij, Starting at `date`"
srun -n 1 ./latlon2ij
echo "Ending at `date`"
