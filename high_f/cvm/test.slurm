#!/bin/bash
#SBATCH -J test_cvm         # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 225              # Total # of nodes
##SBATCH -n 7200           # Total # of mpi tasks
#SBATCH -t 00:05:00        # Run time (hh:mm:ss)
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
srun -n 7200 $UCVMBIN/ucvm2mesh_mpi -f ./test.conf
#srun -n 1620 --exclusive $UCVMBIN/ucvm2mesh_mpi_layer -f ./ucvm_la_habra_small_8m_3456.conf -l 1 -c 5 &
#srun -n 1620 --exclusive $UCVMBIN/ucvm2mesh_mpi_layer -f ./ucvm_la_habra_small_8m_3456.conf -l 6 -c 5 &
#srun -n 1620 --exclusive $UCVMBIN/ucvm2mesh_mpi_layer -f ./ucvm_la_habra_small_8m_3456.conf -l 11 -c 5 &
#srun -n 1620 --exclusive $UCVMBIN/ucvm2mesh_mpi_layer -f ./ucvm_la_habra_small_8m_3456.conf -l 16 -c 5 &
#srun -n 1620 --exclusive $UCVMBIN/ucvm2mesh_mpi_layer -f ./ucvm_la_habra_small_8m_3456.conf -l 21 -c 5 &
#srun -n 1620 --exclusive $UCVMBIN/ucvm2mesh_mpi_layer -f ./ucvm_la_habra_small_8m_3456.conf -l 26 -c 5 &
#wait
echo "ENDING `date`"
