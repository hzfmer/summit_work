#!/bin/bash -l
FIRST=$(qsub dyn.pbs)
echo "$FIRST dyn.pbs"
SECOND=$(qsub -W depend=afterok:$FIRST gmrot.pbs)
echo "$SECOND gmrot.pbs"
THIRD=$(qsub -W depend=afterok:$FIRST pgv2.pbs)
echo "$THIRD pgv2.pbs"
FOURTH=$(qsub -W depend=afterok:$FIRST psr.pbs)
echo "$FOURTH psr.pbs"
FIFTH=$(qsub -W depend=afterok:$FIRST exstat.pbs)
echo "$FIFTH exstat.pbs"
