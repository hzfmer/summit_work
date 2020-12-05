#!/bin/bash
#SBATCH -J test_cvm         # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 72              # Total # of nodes
#SBATCH -t 03:00:00        # Run time (hh:mm:ss)
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)

# Any other commands must follow all #SBATCH directives...
#module list
#pwd
#module load mvapich2
# Launch MPI code...

#export I_MPI_SHM_LMT=shm
#export I_MPI_OFA_DYNAMIC_QPS=0


echo "STARTING `date`"
echo $0
module swap intel gcc
DIR=/ccs/home/hzfmer/scratch/ucvm-19.4.0
UCVMBIN=$DIR/bin
srun -n 2592 $UCVMBIN/ucvm2mesh_mpi -f ./ucvm_la_habra_large_0.conf;
echo "ENDING `date`"
