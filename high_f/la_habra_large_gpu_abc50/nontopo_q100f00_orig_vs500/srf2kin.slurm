#!/bin/bash
#SBATCH -J srf2kin       # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 4              # Total # of nodes
##SBATCH -n 125           # Total # of mpi tasks
#SBATCH -t 00:02:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

echo $0
echo "Running srf2kin, Starting at `date`"
srun -n 125 ~/srf2kin_gp/srf2kin_andes --nx 125 --nz 125 --nio 5 --dt1 0.001 --dt2 0.001 --nt1 5000 --nt2 5000
echo "Ending at `date`"
