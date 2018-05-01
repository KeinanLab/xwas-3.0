#!/bin/sh

# runs qc for X chromosome after imputation

paramfile=$1

## Program will now parse relevant arguments from the param file ###
#------------------------------------------------------------------#
echo "1: Parsing paremeter file"

fname=$( awk '$1=="filename" {print $2}' $paramfile ) # name of dataset file (without the extension)
xwasloc=$( awk '$1=="xwasloc" {print $2}' $paramfile ) # location of extended plink executable
exclind=$( awk '$1=="exclind" {print $2}' $paramfile ) # 0 if no file exists, 1 if file exists
bld=$( awk '$1=="build" {print $2}' $paramfile ) # build of the dataset (18-19)
alpha=$( awk '$1=="alpha" {print $2}' $paramfile ) # alpha threshold to set
binaryplink=$( awk '$1=="plinkformat" {print $2}' $paramfile ) # ped = ped/map format, bed = binary plink format
minmaf=$( awk '$1=="maf" {print $2}' $paramfile ) # remove variants with MAF < minmaf
mindthresh=$( awk '$1=="missindiv" {print $2}' $paramfile ) # filter for missingness per individual
genothresh=$( awk '$1=="missgeno" {print $2}' $paramfile ) # filter for missingness per genotype
quant=$( awk '$1=="quant" {print $2}' $paramfile ) #is this a quantitative trait (either 0 or 1)

if [ $exclind -eq 0 ]; then
    echo -n "" > ${fname}_exclind.remove
    cp ${fname}_exclind.remove ${fname}_male_exclind.remove
    cp ${fname}_exclind.remove ${fname}_female_exclind.remove
fi

if [ $binaryplink = "ped" ] ; then
    ${xwasloc}/xwas --file ${fname} --make-bed --out ${fname} --silent
fi

totsnps=$(wc -l ${fname}.bim | awk '{print $1}')
bonf=$(echo "scale=20;$alpha/$totsnps" | bc)

# Setting PAR locations
if [ $bld -eq 19 ]; then
    echo  '23 60001 2699520 par1\n' > pars.txt
    echo  '23 154931044 155260560 par2' >> pars.txt
fi

if [ $bld -eq 18 ]; then
    echo  '23 1 2709520 par1\n' > pars.txt
    echo  '23 154584238 154913754 par2' >> pars.txt
fi

if [ $bld -eq 17 ]; then
    echo  '23 1 2692881 par1\n' > pars.txt
    echo  '23 154494747 154824264 par2' >> pars.txt
fi

#------------------- Quality Controlling ---------------------------#
#-------------------------------------------------------------------#

# remove pseudo-autosomal regions
echo "2: Removing pseudo-autosomal regions"
${xwasloc}/xwas --bfile ${fname} --make-bed --out ${fname}_psar0 --exclude pars.txt --range --silent

# check for wrong sex individuals (do automatically if > 1000 snps on the X chromosome) Otherwise do this manually
echo "3: Checking wrong-sex individuals"
xtot=`awk '$1==23 {print $0}' ./${fname}_psar0.bim | wc -l`
if [ $xtot -gt 999 ]; then
    ${xwasloc}/xwas --bfile ${fname}_psar0 --check-sex --out ${fname}_psar_sexcheck --silent
    cat ${fname}_psar_sexcheck.sexcheck | grep PROBLEM | awk '{print $1,$2,"sex"}' > ${fname}_sex_exclind.remove
fi

${xwasloc}/xwas --bfile ${fname}_psar0 --make-bed --out ${fname}_psar --remove ${fname}_sex_exclind.remove --silent
malenum=`awk   '$5 == 1' ${fname}_psar.fam | wc -l`;
femalenum=`awk '$5 == 2' ${fname}_psar.fam | wc -l`;

if [ $malenum -eq 0 ] || [ $femalenum -eq 0 ]; then
    echo "Only one gender exist. Male # is $malenum; Female # is $femalenum"
    echo "QC is only performed in one gender"
    if [ $quant -eq 1 ]; then
	echo "4: Quality Control for Quantitative Traits"
	echo "4.1: HWE"
	${xwasloc}/xwas --bfile ${fname}_psar --hardy --out ${fname}_1gender_hwe --silent
	cat ${fname}_1gender_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_1gender_snp.exclude
	echo "4.2: MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_psar --make-bed --out ${fname}_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_1gender_snp.exclude --silent
    else
	echo "4: Quality Control for Binary Traits"
	echo "4.1: HWE"
	${xwasloc}/xwas --bfile ${fname}_psar --hardy --out ${fname}_1gender_hwe --filter-controls --silent
	cat ${fname}_1gender_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_1gender_snp.exclude

	echo "4.2: Correlation between missingness and phenotype"
	${xwasloc}/xwas --bfile ${fname}_psar --test-missing --out ${fname}_1gender_mcc --silent
	awk -v bf=$bonf '$5<bf {print $2}' ${fname}_1gender_mcc.missing >> ${fname}_1gender_snp.exclude
	echo "4.3: MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_psar --make-bed --out ${fname}_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_1gender_snp.exclude --silent
    fi
else
    echo "QC is performed separately for each gender"
    if [ $quant -eq 1 ]; then
	echo "4: Quality Control for Quantitative Traits"
	echo "4.1: Quality Control for Male"
	${xwasloc}/xwas --bfile ${fname}_psar --filter-males --make-bed --out ${fname}_male --silent
	
	echo "4.1.1: HWE"
	${xwasloc}/xwas --bfile ${fname}_male --hardy --out ${fname}_male_hwe --silent 
	cat ${fname}_male_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_male_snp.exclude

	echo "4.1.2: MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_male --make-bed --out ${fname}_male_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_male_snp.exclude --silent

	echo "4.2: Quality Control for Female"
	${xwasloc}/xwas --bfile ${fname}_psar --filter-females --make-bed --out ${fname}_female --silent
	
	echo "4.2.1: HWE"
	${xwasloc}/xwas --bfile ${fname}_female --hardy --out ${fname}_female_hwe --silent 
	cat ${fname}_female_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_female_snp.exclude

	echo "4.2.2: MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_female --make-bed --out ${fname}_female_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_female_snp.exclude --silent

	echo "5: Merging male and female QC"
	perl ${xwasloc}/compare_SNPs_in_male_and_female.pl  ${fname}_female_qc1.bim  ${fname}_male_qc1.bim  ${fname}_sex_diff_snps.diff
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_female_qc1_tmp --silent
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_male_qc1_tmp --silent
	${xwasloc}/xwas --bfile ${fname}_female_qc1_tmp --bmerge ${fname}_male_qc1_tmp.bed ${fname}_male_qc1_tmp.bim ${fname}_male_qc1_tmp.fam --make-bed --out ${fname}_qc1 --silent
    else 
	echo "4: Quality Control for Binary Traits"
	echo "4.1: Quality Control for Male"

	${xwasloc}/xwas --bfile ${fname}_psar --filter-males --make-bed --out ${fname}_male --silent
	
	echo "4.1.1: HWE"
	${xwasloc}/xwas --bfile ${fname}_male --hardy --out ${fname}_male_hwe --filter-controls --silent 
	cat ${fname}_male_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_male_snp.exclude
	
	echo "4.1.2: Correlation between missingness and phenotype"
	${xwasloc}/xwas --bfile ${fname}_male --test-missing --out ${fname}_male_mcc --silent  
	awk -v bf=$bonf '$5<bf {print $2}' ${fname}_male_mcc.missing >> ${fname}_male_snp.exclude

	echo "4.1.3: MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_male --make-bed --out ${fname}_male_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_male_snp.exclude --silent

	echo "4.2: Quality Control for Female"
	${xwasloc}/xwas --bfile ${fname}_psar --filter-females --make-bed --out ${fname}_female --silent
	
	echo "4.2.1: HWE"
	${xwasloc}/xwas --bfile ${fname}_female --hardy --out ${fname}_female_hwe --filter-controls --silent 
	cat ${fname}_female_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_female_snp.exclude
	
	echo "4.2.2: Correlation between missingness and phenotype"
	${xwasloc}/xwas --bfile ${fname}_female --test-missing --out ${fname}_female_mcc --silent  
	awk -v bf=$bonf '$5<bf {print $2}' ${fname}_female_mcc.missing >> ${fname}_female_snp.exclude

	echo "4.2.3: MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_female --make-bed --out ${fname}_female_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_female_snp.exclude --silent
	
	echo "5: Merging male and female QC"
	perl  ${xwasloc}/compare_SNPs_in_male_and_female.pl  ${fname}_female_qc1.bim  ${fname}_male_qc1.bim  ${fname}_sex_diff_snps.diff
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_female_qc1_tmp --silent
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_male_qc1_tmp --silent
	${xwasloc}/xwas --bfile ${fname}_female_qc1_tmp --bmerge ${fname}_male_qc1_tmp.bed ${fname}_male_qc1_tmp.bim ${fname}_male_qc1_tmp.fam --make-bed --out ${fname}_qc1 --silent
    fi
fi

# Significant diff in MAF between males and females only for qualitative traits and only in controls
if [ $quant -eq 1 ]; then
    echo "6: QC on female-male merged dataset: MAF; Individual Missingness; Genotyp Missingness"
    ${xwasloc}/xwas --bfile ${fname}_qc1  --make-bed --out ${fname}_final_x  --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh}  --silent
else
    echo "6: QC on female-male merged dataset: MAF; Individual Missingness; Genotyp Missingness"
    ${xwasloc}/xwas --bfile ${fname}_qc1  --make-bed --out ${fname}_final  --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh}  --silent

    echo "7: Significant difference in MAF between males and females (Only for binary traits and Only in Controls)"
    totx=$(awk '$1=="X" || $1==23 {print $0}' ${fname}_final.bim | wc -l)
    bonfx=$(echo "scale=20;$alpha/$totx" | bc)
    ${xwasloc}/xwas --bfile ${fname}_final  --xwas --make-bed --out ${fname}_final_x --freqdiff-x ${bonfx} --silent
fi

# Clean up
rm ${fname}_female* ${fname}_male* ${fname}_psar* ${fname}_qc1* ${fname}_sex*
rm pars.txt *.exclude *.remove
