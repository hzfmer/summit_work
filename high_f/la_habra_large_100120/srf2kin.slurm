#!/bin/bash
#SBATCH -J srf2kin       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 4              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 00:02:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo $0
echo "Running srf2kin, Starting at `date`"
srun -n 125 ~/srf2kin_gp/srf2kin_andes 125 125 5 0.001 0.0004 5000 12500
echo "Ending at `date`"
