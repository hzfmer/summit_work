#!/bin/bash -l

while true
do
  sleep 10m
  if [ -f sliprate.bin -a ! -f sliprate_1hz.bin ]; then
    break
  fi

ZERO=$(qsub srf2kin.pbs)
echo "$ZERO srf2kin.pbs"
./pre-run-srcpart
FIRST=$(qsub -W depend=afterok:$ZERO srcpart.pbs)
echo "$FIRST srcpart.pbs"
SECOND=$(qsub -W depend=afterok:$FIRST dyn.pbs)
echo "$SECOND dyn.pbs"
THIRD=$(qsub -W depend=afterok:$SECOND pgv2.pbs)
echo "$THIRD pgv2.pbs"
FIFTH=$(qsub -W depend=afterok:$SECOND extrts.pbs)
echo "$FIFTH exstat.pbs"
