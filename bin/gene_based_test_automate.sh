#!/bin/bash

# implements truncated tail method for gene-based GWAS, accounting for LD

paramfile=$1
## Program will now parse relevant arguments from the param file ###
#------------------------------------------------------------------#
echo "Parsing arguments"

#run script with parameters
fname=$( awk '$1=="filename" {print $2}' $paramfile ) # name of dataset file (without the extension)
xwasloc=$( awk '$1=="xwasloc" {print $2}' $paramfile ) # location of xwas executable
genetestloc=$( awk '$1=="genescriptloc" {print $2}' $paramfile ) # location of the gene-based R script.
genelist=$( awk '$1=="genelistname" {print $2}' $paramfile ) # list of genes (with location)
pvfolder=$( awk '$1=="pvfolder" {print $2}' $paramfile ) # file with snp names [column1] and p-values [column2]
buff=$( awk '$1=="buffer" {print $2}' $paramfile ) # buffer around a gene (in base pairs)
#------------------------------------------------------------------#

chr="23"
pvcol=2
snpcol=1

# make a bfile if it doesn't exist
if [ ! -f ${fname}.bed ];
then
  ${xwasloc}xwas --file ${fname} --make-bed --out ${fname} --noweb --silent
fi

# for each *_temp file in the pvfolder
for f in ${pvfolder}*_temp.*;
do

	# create a gene-based p-value text file
	outputname=${f/_temp/_gene_results} #create output file by removing "_temp*" part and add "_gene_results.txt"
	echo "Writing results to [ ${outputname} ]"
	echo "gene reps tail_p_value prod_p_value" > ${outputname}

	# because plink doesn't process bp past the last snp...
	firstbp=$(head -n 1 ${fname}.bim | awk '{print $4}')
	lastbp=$(tail -1 ${fname}.bim | awk '{print $4}')

	while read -r chr stbp endbp gene
	do

  	        echo "Testing $gene"

		# buffers
		buffstart=$(echo "$stbp-$buff" | bc)
		buffstop=$(echo "$endbp+$buff" | bc)

		# make sure regions don't overshoot the boundaries...
		if [ $buffstart -lt $firstbp ]; then
			buffstart=$firstbp
		fi
		if [ $buffstop -gt $lastbp ]; then
			buffstop=$lastbp
		fi

		# only do this if the gene is within the region of our data
		if [ $buffstart -lt $lastbp ] && [ $buffstop -gt $firstbp ]; then
			# making the smaller plink file
			# get a list of snps in this gene
			awk '($1 == '"$chr"') && ($4 > '"$buffstart"') && ($4 < '"$buffstop"') {print $2}' ${fname}.bim > ${fname}_${gene}.snps

			numsnps=$(wc -l ${fname}_${gene}.snps | awk '{print $1}' )
			if [ $numsnps -ne 0 ]; then

				# extract the snps
				${xwasloc}/xwas --bfile ${fname} --extract ${fname}_${gene}.snps --out temp --make-bed --noweb --silent

				# limit to only SNPs within buffer region of the gene
				awk 'FNR==NR {f1[$2]; next} ($'${snpcol}' in f1) {print $'${snpcol}',$'${pvcol}'}' temp.bim ${f} > temppv.txt

				# temp snp list
				awk '{print $1}' temppv.txt > temp.snp

				# get the LD matrix
				${xwasloc}/xwas --bfile temp --matrix --r --extract temp.snp --noweb --out ${fname}_ld --silent

				# replacing nan values with 0
			        perl -i -pe 's/nan/0/g' ${fname}_ld.ld

				# do the trunctated test
	                        Rscript ${genetestloc}/truncgene.R ${gene} ${outputname} ${fname}

	                        rm ${fname}_ld.* temp.*
			fi
			rm ${fname}_${gene}.snps
		fi
	done <"$genelist"

	sort -nk 3 ${outputname} > ${outputname}.sort
	rm temppv.txt ${outputname}
        mv ${outputname}.sort ${outputname}
done
