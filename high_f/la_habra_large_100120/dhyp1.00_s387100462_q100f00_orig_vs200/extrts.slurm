#!/bin/bash
#SBATCH -J extrts       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 1              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 01:40:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

# Require IN3D.out
echo "Running extrts, Starting at `date`"
srun -n 1  ~/file_back/programs/awp_tools/bin_andes/extrts
wait
module load python/3.7.0-anaconda3-2018.12
# python pickle_sites.py ../la_habra_small_statlist_3456.idx;
echo "Ending at `date`"
