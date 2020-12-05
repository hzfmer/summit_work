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
COUNT=$(($(ls ./S?_0_00*0 | wc -l) / 3))
STEP=$((NT / COUNT))
echo "COUNT=$COUNT, STEP=$STEP for each direction"
echo S${F}_0_$(printf %07d $((i * STEP)))

for F in X Y Z; do
  FILE=S${F}_0_
  echo ${FILE} && touch ${FILE}
  [ -f ${FILE} ] && rm ${FILE}
  for i in {1..$COUNT}; do
    cat S${F}_0_$(printf %07d $((i * STEP))) >> ${FILE}
  done
done
cd ..

