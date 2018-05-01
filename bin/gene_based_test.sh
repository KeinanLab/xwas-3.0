#!/bin/sh

# implements truncated tail method for gene-based GWAS, accounting for LD

paramfile=$1

## Program will now parse relevant arguments from the param file ###
#------------------------------------------------------------------#
echo "Parsing arguments"

fname=$( awk '$1=="filename" {print $2}' $paramfile ) # name of dataset file (without the extension)
xwasloc=$( awk '$1=="xwasloc" {print $2}' $paramfile ) # location of xwas executable
genetestloc=$( awk '$1=="genescriptloc" {print $2}' $paramfile ) # location of the gene-based R script. 
genelist=$( awk '$1=="genelistname" {print $2}' $paramfile ) # list of genes (with location)
pvlist=$( awk '$1=="assocfile" {print $2}' $paramfile ) # file with snp names [column1] and p-values [column2]
buff=$( awk '$1=="buffer" {print $2}' $paramfile ) # buffer around a gene (in base pairs)
outputname=$( awk '$1=="output" {print $2}' $paramfile ) # name of output file
numindiv=$( awk '$1=="numindiv" {print $2}' $paramfile ) # number of individuals to estimate LD
#------------------------------------------------------------------#

chr="23"
pvcol=2
snpcol=1

# first creat a custom genotype file made of random controls. 
awk '$6==1 {print $0}' ${fname}.fam | head -${numindiv} | awk '{print $1"\t"$2}' > ${fname}_rand200controls.keep
# make a bfile
${xwasloc}/xwas --bfile ${fname} --make-bed --keep ${fname}_rand200controls.keep --out ${fname}_rand200controls --noweb  --silent

# create a gene-based p-value text file
echo "gene reps tail_p_value prod_p_value" > ${outputname}

# because plink doesn't process bp past the last snp...
firstbp=$(head -n 1 ${fname}_rand200controls.bim | awk '{print $4}')
lastbp=$(tail -1 ${fname}_rand200controls.bim | awk '{print $4}')


echo "Testing genes"
while read -r chr stbp endbp gene
do

	echo "Testing : $gene"
	
	# buffers
	buffstart=$(echo "$stbp-$buff" | bc)
	buffstop=$(echo "$endbp+$buff" | bc)
	# make sure regions don't overshoot the boundaries...
	if [ $buffstart -lt $firstbp ]
	then
		buffstart=$firstbp
	fi
	if [ $buffstop -gt $lastbp ]
	then
		buffstop=$lastbp
	fi
	
	# only do this if the gene is within the region of our data
	if [ $buffstart -lt $lastbp ] && [ $buffstop -gt $firstbp ]
	then
		# making the smaller plink file
		# get a list of snps in this gene
		awk '($1 == '"$chr"') && ($4 > '"$buffstart"') && ($4 < '"$buffstop"') {print $2}' ${fname}_rand200controls.bim > ${fname}_${gene}.snps	

		numsnps=$(wc -l ${fname}_${gene}.snps | awk '{print $1}' )
		if [ $numsnps -gt 0 ]
		then

			# extract the snps
			${xwasloc}/xwas --bfile ${fname}_rand200controls --extract ${fname}_${gene}.snps --out temp --make-bed --silent --noweb

			# limit to only SNPs within buffer region of the gene
			awk 'FNR==NR {f1[$2]; next} ($'${snpcol}' in f1) {print $'${snpcol}',$'${pvcol}'}' temp.bim ${pvlist} > temppv.txt

			# temp snp list
			awk '{print $1}' temppv.txt > temp.snp
		
			# get the LD matrix
			${xwasloc}/xwas --bfile temp --matrix --r --extract temp.snp --noweb --out ${fname}_ld --silent
			# replacing nan values with 0
			perl -i -pe 's/nan/0/g' ${fname}_ld.ld
		
			# do the trunctated test
                        Rscript ${genetestloc}/truncgene.R ${gene} ${outputname} ${fname}

                        rm ${fname}_ld.ld ${fname}_ld.log temp.*
		fi
		rm ${fname}_${gene}.snps
	fi
done <"$genelist"

sort -nk 3 ${outputname} > ${outputname}.sort

rm temppv.txt ${fname}_rand200controls.*

