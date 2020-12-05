#!/bin/bash
#SBATCH -J test_cvm         # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 1              # Total # of nodes (16 / 32)
#SBATCH -t 00:03:00        # Run time (hh:mm:ss)
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

# Any other commands must follow all #SBATCH directives...

echo "STARTING `date`"
echo $0
DIR=/ccs/home/nanwang/ucvm-19.4.0
UCVMBIN=$DIR/bin
srun -n 16 $UCVMBIN/ucvm2mesh_mpi -f test.conf
echo "ENDING `date`"
