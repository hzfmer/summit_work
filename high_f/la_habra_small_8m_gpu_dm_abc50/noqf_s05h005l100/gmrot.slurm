#!/bin/bash
#SBATCH -J gmrot       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 18              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 00:40:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

## 576 divides ny
echo "Running gmrot, Starting at `date`"
#srun -n 576  ~/file_back/programs/awp_tools/bin_rhea/gmrot_mpi 1.0 2.0 3.0 4.0 5.0
srun -n 576  ~/file_back/programs/awp_tools/bin_rhea/gmrot_mpi 0.5 1 2.5 3.5
echo "Ending at `date`"
