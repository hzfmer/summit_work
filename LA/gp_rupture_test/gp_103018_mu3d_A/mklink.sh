#!/bin/bash -l

DIR=output_sfc/

nt=$1
echo ${nt}
cd $DIR
for file in `find -type f`; do
  ln -s ${file} ${file}_$(printf %07d $nt)
done

