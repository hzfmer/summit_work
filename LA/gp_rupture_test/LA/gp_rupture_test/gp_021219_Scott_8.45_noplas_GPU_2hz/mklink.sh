#!/bin/zsh
set -e

INPUT=param.sh
TMAX=$(awk -F 'TMAX' '{print $2}' ${INPUT} | awk '{print $1}')
DT=$(awk -F 'DT' '{print $2}' ${INPUT} | awk '{print $1}')
echo ${TMAX} ${DT}
NT=$(printf %.0f $((TMAX / DT)))

echo "Linking SX/SY/SZ, NT=${NT}"
number=$(printf %07d ${NT})

cd output_sfc
for F in X Y Z; do
  if [ -f S${F}${number} ]; then
    rm S${F}${number}
  fi
  ln -s S$F S${F}${number}
  echo "Linking S$F"
done
cd ..


