#!/bin/bash -l

rm sliprate.bin
rm sliprate_2hz.bin
rm momrate.dat

lfs setstripe -c 12 sliprate.bin
lfs setstripe -c 12 sliprate_2hz.bin
lfs setstripe -c 48 momrate.dat
lfs setstripe -c 48 output_sfc/SX
lfs setstripe -c 48 output_sfc/SY
lfs setstripe -c 48 output_sfc/SZ

echo "Running: python comp_sliprate.py"
python comp_sliprate.py &
while true
do 
  sleep 5m
  if [ -f sliprate.bin ]; then
    A=$(stat --printf="%s" sliprate.bin)
    if [ $A == "18008320000" ]; then
      break
    fi
  fi
done

python filter_sr.py &
while true
do 
  sleep 5m
  if [ -f sliprate_2hz.bin ]; then
    A=$(stat --printf="%s" sliprate_2hz.bin)
    if [ $A == "18008320000" ]; then
      break
    fi
  fi
done
mv sliprate_2hz.bin sliprate.bin

while true
do
  sleep 1m
  if [ -f sliprate.bin -a ! -f sliprate_2hz.bin ]; then
    break
  fi
done

ZERO=$(qsub srf2kin.pbs)
echo "$ZERO srf2kin.pbs"
./pre-run-srcpart
FIRST=$(qsub -W depend=afterok:$ZERO srcpart.pbs)
echo "$FIRST srcpart.pbs"
SECOND=$(qsub -W depend=afterok:$FIRST dyn.pbs)
echo "$SECOND dyn.pbs"
#THIRD=$(qsub -W depend=afterok:$SECOND pgv2.pbs)
#echo "$THIRD pgv2.pbs"
FIFTH=$(qsub -W depend=afterok:$SECOND extrts.pbs)
echo "$FIFTH extrts.pbs"
SIXTH=$(qsub -W depend=afterok:$SECOND gmrot.pbs)
echo "$SIXTH gmrot.pbs"
