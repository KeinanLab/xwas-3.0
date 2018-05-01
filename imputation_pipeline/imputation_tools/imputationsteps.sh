#!/bin/sh

fname=$1
wname=$2
build=$3

# convert to ped/map format
./xwas --noweb --bfile ${fname} --recode --out $wname

# liftover 
if [ ${build} != 19 ]
then
    echo "Genome build other than 19 detected, must use liftOver"
    ./liftme_fast.sh ${wname}.map ${wname}.hg19.map ${build}
else
    cp ${wname}.map  ${wname}.hg19.map
fi

# convert to impute2 format
echo "Converting to impute2 format..."
./gtool -P --ped ${wname}.ped --map ${wname}.hg19.map
