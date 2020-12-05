#!/bin/bash -l

ZERO=$(qsub srf2kin.pbs)
echo "$ZERO srf2kin.pbs"
./pre-run-srcpart
FIRST=$(qsub -W depend=afterok:$ZERO srcpart.pbs)
echo "$FIRST srcpart.pbs"
SECOND=$(qsub -W depend=afterok:$FIRST dyn.pbs)
echo "$SECOND dyn.pbs"
THIRD=$(qsub -W depend=afterok:$SECOND pgv2.pbs)
echo "$THIRD pgv2.pbs"
FIFTH=$(qsub -W depend=afterok:$SECOND exstat.pbs)
echo "$FIFTH exstat.pbs"
