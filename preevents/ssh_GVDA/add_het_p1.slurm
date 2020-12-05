#!/bin/bash
#SBATCH -J add_het_p0  # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
##SBATCH -N 47              # Total # of nodes
##SBATCH -n 1500           # Total # of mpi tasks
#SBATCH -N 15              # Total # of nodes
##SBATCH -n 464           # Total # of mpi tasks
#SBATCH -t 00:05:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo "Running add_het_p1, Starting at `date`"
srun -n 464 ./add_het_p1
#srun -n 1500 ./add_het_p1
echo "Ending at `date`"
