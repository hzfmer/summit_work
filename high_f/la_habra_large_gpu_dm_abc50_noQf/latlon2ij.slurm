#!/bin/bash
#SBATCH -J latlon2ij       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 1              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 00:35:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo "Running latlon2ij_2D, Starting at `date`"
srun -n 1 ./latlon2ij_2D 9504 7020 surf.grid la_habra_large_statlist.txt test_latlon2D.txt
echo "Ending at `date`"
