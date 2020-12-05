#!/bin/bash
# Begin LSF Directives
#BSUB -J high_f_sm_1layer
#BSUB -P geo112
#BSUB -W 01:20
#BSUB -nnodes 96
#BSUB -alloc_flags "gpumps"
##BSUB -alloc_flags maximizegpfs
#BSUB -o run_%J.out
#BSUB -e run_%J.err
##BSUB -N

export OMP_NUM_THREADS=1
cd $LS_SUBCWD

module unload darshan-runtime
module load cuda

cat $0
args=$(cat param_origsrc.sh)
echo $args

#jsrun -n 72 -a 3 -c 3 -g 1 -r 6 -d cyclic /ccs/home/hzfmer/scratch/AWP-ODC-YFSWP/pmcl3d $args
jsrun -n 576 -a 3 -c 3 -g 1 -r 6 -d cyclic /ccs/home/hzfmer/awp_highf/pmcl3d_qs100_vs200 $args
