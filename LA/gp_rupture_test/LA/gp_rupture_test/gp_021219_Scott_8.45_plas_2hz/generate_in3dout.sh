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
cp IN3D IN3D.out
cp IN3D IN3D_gmrot.out
sed -i "s/^.*WRITE_STEP[^0-9]/${NT} WRITE_STEP/" IN3D_gmrot.out
