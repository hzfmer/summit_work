#!/bin/zsh

#CASE=$(echo ${PWD##*/} | cut -d '_' -f 3-)
CASE=$(echo ${PWD##*/})
echo $CASE
if [ -f param.sh ]; then
  INPUT=param.sh
else
  INPUT=run.lsf
fi
TMAX=$(awk -F 'TMAX' '{print $2}' ${INPUT} | awk '{print $1}')
DT=$(awk -F 'DT' '{print $2}' ${INPUT} | awk '{print $1}')
NT=$(printf %.0f $((TMAX / DT)))

cd output_sfc
COUNT=$(($(ls ./S* | wc -l) / 3))
STEP=$((NT / COUNT))
echo "COUNT=$COUNT, STEP=$STEP for each direction"

for F in X Y Z; do
  FILE=v$(tr '[A-Z]' '[a-z]' <<< $F)_${CASE}.bin
  echo $FILE
  [ -f $FILE ] && rm $FILE
  for i in {1..$COUNT}; do
    cat S${F}_0_$(printf %07d $((i * STEP))) >> $FILE
  done
  tomac.sh -d /Users/zhh076/work/preevents/data $FILE
done
cd ..

