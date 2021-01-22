#!/bin/bash
#SBATCH -J ssh_8m_h005l100  # Job name
#SBATCH -o %j.out          # Name of stdout output file
#SBATCH -e %j.err          # Name of stderr error file
#SBATCH -N 1              # Total # of nodes
##SBATCH -n 1           # Total # of mpi tasks
#SBATCH -t 4:00:00        # Run time (hh:mm:ss)
#SBATCH -p gpu
##SBATCH --mail-type=fail    # Send email at begin and end of job
##SBATCH --mail-user=zhh076@ucsd.edu
#SBATCH -A geo112       # Project/Allocation name (req'd if you have more than 1)
cd $SLURM_SUBMIT_DIR
module load gcc fftw/3.3.8

echo $PWD
echo $SLURM_SUBMIT_DIR
echo "Running mkssh, Starting at `date`"
echo $0


#srun -n 1 ./ssh_generate -m ssh_20m_large_single_r05h005l100.out --d1 20 --hurst 0.05 --l1 100 --st23 5 --n1 500 --n2 9504 --n3 7020
srun -n 1 ./ssh_generate -m ssh_20m_large_single_r05h005l500.out --d1 20 --hurst 0.05 --l1 500 --st23 5 --n1 500 --n2 9504 --n3 7020
#mpirun -n 1 ssh_generate_dm -m ssh_20m_dm.out --d1 20 --hurst 0.1 --l1 100 --st23 5 --n1 500 --n2 9000 --n3 6750 --dmi 360
#mpirun -n 1 ssh_generate_dm -m ssh_20m_sm.out --d1 20 --hurst 0.1 --l1 100 --st23 5 --n1 500 --n2 1400 --n3 1400 --dmi 360
#srun -n 1 ssh_generate_dm -m ssh_8m_r05h010l100.out --d1 8 --hurst 0.10 --l1 100 --st23 5 --n1 1250 --n2 3456 --n3 3456 --dmi 108
#srun -n 1 ssh_generate_dm -m ssh_8m_r05h015l100.out --d1 8 --hurst 0.15 --l1 100 --st23 5 --n1 1250 --n2 3456 --n3 3456 --dmi 108
#srun -n 1 ssh_generate_dm -m ssh_8m_r05h005l200.out --d1 8 --hurst 0.05 --l1 200 --st23 5 --n1 1250 --n2 3456 --n3 3456 --dmi 108
#srun -n 1 ssh_generate_dm -m ssh_8m_r05h005l500.out --d1 8 --hurst 0.05 --l1 500 --st23 5 --n1 1250 --n2 3456 --n3 3456 --dmi 108
#srun -n 1 ssh_generate_dm -m ssh_8m_r05h010l500.out --d1 8 --hurst 0.10 --l1 500 --st23 5 --n1 1250 --n2 3456 --n3 3456 --dmi 108
#srun -n 1 ssh_generate_dm -m ssh_8m_r05h005l050.out --d1 8 --hurst 0.05 --l1 50 --st23 5 --n1 1250 --n2 3456 --n3 3456 --dmi 108
#srun -n 1 ssh_generate_dm -m ssh_8m_r05h005l500.out --d1 8 --hurst 0.05 --l1 500 --st23 5 --n1 1250 --n2 3456 --n3 3456 --dmi 108
#srun -n 1 ssh_generate_dm -m ssh_8m_large_r05h005l100.out --d1 8 --hurst 0.05 --l1 100 --st23 5 --n1 1250 --n2 19440 --n3 14904 --dmi 184
#srun -n 1 ssh_generate_dm -m ssh_8m_large_r05h005l500.out --d1 8 --hurst 0.05 --l1 500 --st23 5 --n1 1250 --n2 19440 --n3 14904 --dmi 184
echo "Ending at `date`"
