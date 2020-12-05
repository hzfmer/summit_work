#!/bin/bash
#SBATCH -J add_het_p1  # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 80              # Total # of nodes
#SBATCH -n 1250            # Total # of mpi tasks
#SBATCH -t 00:35:00        # Run time (hh:mm:ss)
#SBATCH --mem=0            # Use all available memory 
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)
module unload Darshan
echo "Running add_het_p1, Starting at `date`"
# -n = nz of the hom grid
#srun -n 384 ./add_het_p1
# small
# srun -n 1400 ./add_het_p1
# ext-large
# srun -n 3072 ./add_het_p1
# large-8m
srun -n 1250 -N 80 ./add_het_p1_simple
echo "Ending at `date`"
