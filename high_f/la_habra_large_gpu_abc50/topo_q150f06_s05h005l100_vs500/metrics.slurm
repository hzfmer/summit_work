#!/bin/bash
#SBATCH -J metrics       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 30              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 06:00:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo "Running metrics, Starting at `date`"
srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 1 -n 13 
#wait
#srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 4 -n 13
#wait
#srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 2.5 -n 13
#srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 2.5 -h 5.0 -n 13
#wait
#srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 2 -n 13
#wait
#srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 3 -n 13
#wait
#srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -f -l 0.15 -h 5 -n 13
#wait
#srun -n 960  ~/file_back/programs/awp_tools/bin_andes/metrics -n 13
#wait

echo "Ending at `date`"
