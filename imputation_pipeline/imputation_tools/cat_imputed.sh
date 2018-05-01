#!/bin/sh

# cats the imputed files 
fname=$1
numfiles=$(echo $2-1 | bc)
thresh=$3 # threshold, below which snps are removed because they dont' have good quality score [used 0.5]
samp=$4 # sample file w. location
outname=$5
cmd='cat'

if [ -f ${fname}_lowqsnps.txt ]; then
    rm ${fname}_lowqsnps.txt
fi

for i in `seq 0 $numfiles`;
do
	cmd=$(echo $cmd $fname.$i)
	awk -v cutter=$thresh '$5 < cutter {print $2}' ${fname}.${i}_info >> ${fname}_lowqsnps.txt

done

$cmd > ${fname}.imp
./gtool -S --g ${fname}.imp --s ${samp} --exclusion ${fname}_lowqsnps.txt --og ${fname}.imp.qc 
./gtool -G --g ${fname}.imp.qc --s ${samp} --phenotype phenotype --chr 23 --ped ${fname}.imp.qc.ped --map ${fname}.imp.qc.map --snp

# do some final qc (missingness per snp and sig sex differences in maf
alpha=0.05
totsnps=$(wc -l ${fname}.imp.qc.map | awk '{print $1}')
bonf=$(echo "scale=20;$alpha/$totsnps" | bc)

## revised by Kaixiong
./xwas --noweb --file ${fname}.imp.qc --make-bed --out ${outname}  --missing-genotype N
