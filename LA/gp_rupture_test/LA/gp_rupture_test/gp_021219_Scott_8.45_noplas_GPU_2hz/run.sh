#!/bin/bash -l

python generate_in3dout.py @param.sh &

# rm sliprate.bin
# rm sliprate_2hz.bin
rm momrate.dat

#python comp_sliprate.py &
#echo "Running: python comp_sliprate.py"
#while true
#do
#  if [ -f sliprate.bin ]; then
#    A=$(stat --printf="%s" sliprate.bin)
#    echo $A
#    if [ $A == "18008320000" ]; then
#      break
#    fi
#  fi
#  sleep 1m
#done

#python filter_sr.py &
#while true
#do 
#  if [ -f sliprate_2hz.bin ]; then
#    A=$(stat --printf="%s" sliprate_2hz.bin)
#    if [ $A == "18008320000" ]; then
#      break
#    fi
#  fi
#  sleep 1m
#done
#mv sliprate_2hz.bin sliprate.bin

while true
do
  if [ -f sliprate.bin -a ! -f sliprate_2hz.bin ]; then
    break
  fi
  sleep 1m
done

ZERO=$(bsub srf2kin.lsf)
echo "$ZERO srf2kin.lsf"
./pre-run-srcpart
FIRST=$(bsub -w ${ZERO//[^0-9]/} srcpart.lsf)
echo "$FIRST srcpart.lsf"
SECOND=$(bsub -w ${FIRST//[^0-9]/} run.lsf)
echo "$SECOND run.lsf"
while true
do
  if [[ $(bjobs -d ${SECOND//[^0-9]/}) == *"DONE"* ]]; then
    ./accum_v.sh
  fi
  sleep 10m
done
THIRD=$(bsub -w ${SECOND//[^0-9]/} pgv2.lsf)
echo "$THIRD pgv2.lsf"
FIFTH=$(bsub -w ${SECOND//[^0-9]/} extrts.lsf)
echo "$FIFTH extrts.lsf"
SIXTH=$(bsub -w ${SECOND//[^0-9]/} gmrot.lsf)
echo "$SIXTH gmrot.lsf"
