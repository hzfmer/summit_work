#!/bin/bash -l

lfs setstripe -c 48 output_sfc/SX
lfs setstripe -c 48 output_sfc/SY
lfs setstripe -c 48 output_sfc/SZ


SECOND=$(qsub -W depend=afterok:9771507 dyn.pbs)
echo "$SECOND dyn.pbs"
#THIRD=$(qsub -W depend=afterok:$SECOND pgv2.pbs)
#echo "$THIRD pgv2.pbs"
FIFTH=$(qsub -W depend=afterok:$SECOND extrts.pbs)
echo "$FIFTH extrts.pbs"
SIXTH=$(qsub -W depend=afterok:$SECOND gmrot.pbs)
echo "$SIXTH gmrot.pbs"
