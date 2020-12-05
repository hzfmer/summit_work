#!/bin/zsh
zmodload zsh/mathfunc

line=$(grep -n -m 1 "TMAX" IN3D)
TMAX=$(echo $line | awk '{print $2}')

line=$(grep -n -m 1 "DT" IN3D)
DT=$(echo $line | awk '{print $2}')

line=$(grep -ni -m 1 "NTISKP" IN3D)
NTISKP=$(echo $line | awk '{print $2}')

echo ${DT}
echo ${TMAX}
echo ${NTISKP}
NT=$((int(TMAX/DT/NTISKP)))
echo ${NT}

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
