#!/bin/bash -l

for F in `find ./ -name "m[0-9]*"`;
do 
    MAG=$(echo $F | cut -d '_' -f1 | cut -c4-)
    SRC=source_"$MAG"_lvz
    yi=$(grep "yi=" "$SRC"/define_stress_F.py)
    yi=$(echo $yi | cut -d $'\n' -f 1 | cut -d ' ' -f 1)
    echo $yi
    yi=$((yi + 1))
    nx=$(grep "nx=" "$SRC"/define_stress_F.py)
    nx=$(echo $nx | cut -d $'\n' -f 1 | cut -d ' ' -f 1)
    echo $nx
    nx=$((nx - 20))
    zi1=$(grep "zi1=" "$SRC"/define_stress_F.py)
    zi1=${zi1#*=}
    cd $F
    echo "21 $nx
    $yi
    1 $zi1
    10" > fault_output.dat
    sed -i 's/^[ \t]*//' fault_output.dat
    cd ..
done


