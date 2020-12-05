#!/bin/bash 

echo $1

if [[ $# != 1 ]]; then
  echo "Usage: ./mklink.sh NT"
  exit 1
else
  echo "Linking SX/SY/SZ"
fi

number=$(printf %07d $1)
echo "Number of NT: ${number}"

cd output_sfc
for F in X Y Z; do
  if [ -f S${F}${number} ]; then
    rm S${F}${number}
  fi
  ln -s S$F S${F}${number}
  echo "Linking S$F"
done
cd ..


