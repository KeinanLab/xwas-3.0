#!/bin/bash

# implements gene_based gene x gene interaction test

function calc {
    count=$(cat $1 | wc -l) # number of genes in list1
    count2=$(cat $2 | wc -l) # number of genes in list2
    for(( i=1; i<=$count;i++)); do
        str1=$(sed -n "$i"p $1)

        gene1stbp=$(echo $str1 | awk '{print $2}')
        gene1endbp=$(echo $str1 | awk '{print $3}')
        gene1=$(echo $str1 | awk '{print $4}')
        gene1chr=$(echo $str1 | awk '{print $1}')

        buffstart1=$(echo "$gene1stbp-$buff" | bc)
        buffstop1=$(echo "$gene1endbp+$buff" | bc)

        firstbp1=$(awk '$1 == '"$gene1chr"' {print $4}' ${fname}.bim | head -1)
        lastbp1=$(awk '$1 == '"$gene1chr"' {print $4}' ${fname}.bim | tail -1)

        if [[ $buffstart1 -lt $firstbp1 ]]; then
            buffstart1=$firstbp1
        fi

        if [[ $buffstop1 -gt $lastbp1 ]]; then
            buffstop1=$lastbp1
        fi

        if [[ $buffstart1 -lt $lastbp1 ]] && [[ $buffstop1 -gt $firstbp1 ]]; then

            if [[ -z "${genelist2// }" ]]; then
                start=`expr $i + 1`
            else
                start=1
            fi

            for(( j=$start;j<=$count2;j++)); do
                str2=$(sed -n "$j"p $2)
                gene2stbp=$(echo $str2 | awk '{print $2}')
                gene2endbp=$(echo $str2 | awk '{print $3}')
                gene2=$(echo $str2 | awk '{print $4}')
                gene2chr=$(echo $str2 | awk '{print $1}')
    
                buffstart2=$(echo "$gene2stbp-$buff" | bc)
                buffstop2=$(echo "$gene2endbp+$buff" | bc)
    
                # because plink doesn't process bp past the last snp...
                firstbp2=$(awk '$1 == '"$gene2chr"' {print $4}' ${fname}.bim | head -1)
                lastbp2=$(awk '$1 == '"$gene2chr"' {print $4}' ${fname}.bim | tail -1)
    
                # make sure regions don't overshoot the boundaries...
                if [[ $buffstart2 -lt $firstbp2 ]]; then
                    buffstart2=$firstbp2
                fi
    
                if [[ $buffstop2 -gt $lastbp2 ]]; then
                    buffstop2=$lastbp2
                fi
    
                # only do this if both genes are within the region of our data
                if [[ $buffstart2 -lt $lastbp2 ]] && [[ $buffstop2 -gt $firstbp2 ]]; then
                    echo "SET_1" > ${fname}.set
                    awk ' ($1 == '"$gene1chr"') && ($4 >= '"$buffstart1"') && ($4 <= '"$buffstop1"') {print $2}' ${fname}.bim > ${fname}_temp1.snps
                    cat ${fname}_temp1.snps >> ${fname}.set
                    echo "END" >> ${fname}.set
                    echo >> ${fname}.set
                    echo "SET_2" >> ${fname}.set
                    awk ' ($1 == '"$gene2chr"') && ($4 >= '"$buffstart2"') && ($4 <= '"$buffstop2"') {print $2}' ${fname}.bim > ${fname}_temp2.snps
                    cat ${fname}_temp2.snps >> ${fname}.set
                    echo "END" >> ${fname}.set
    
                    cat ${fname}_temp1.snps ${fname}_temp2.snps > ${fname}_temp.snps
    
                    numsnps=$(wc -l ${fname}_temp.snps | awk '{print $1}')
                    if [[ $numsnps -ne 0 ]]; then
    	                echo "Writing SNP-pairwise epistasis results for $gene1 and $gene2 to [ ${gene1}_${gene2}.qt ]"

                        ${xwasloc}/xwas --bfile ${fname} --extract ${fname}_temp.snps --out temp --make-bed --silent --noweb

                        if [[ -z "$covarfile" ]]
                        then
                            ${xwasloc}/xwas --noweb --xwas --bfile temp --xepi --set ${fname}.set --xepi1 1 --xepi2 1 --out ${fname} --silent
                        else
                            ${xwasloc}/xwas --noweb --xwas --bfile temp --xepi --covar ${covarfile} --set ${fname}.set --xepi1 1 --xepi2 1 --out ${fname} --silent
                        fi
    
                        # Use R script to calc the four different pvalues (params: pvalue file, corrlation matrix file, summary output file)
                        Rscript ${geneinterloc}/gene_based_inter.R ${fname}.xepi.qt ${fname}.xepi.qt.corr $gene1 $gene2 ${outputname}
    
                        mv ${fname}.xepi.qt ${gene1}_${gene2}.qt
                        rm ${fname}.xepi.qt.* *.snps *.set temp.*
                    fi
                fi
            done
        fi
    done
}

#------------------------------------------------------------------#
#                             MAIN                                 #
#------------------------------------------------------------------#

paramfile=$1

echo "Parsing arguments"

fname=$( awk '$1=="filename" { print $2 }' $paramfile ) # name of the dataset file (without the extension)
xwasloc=$( awk '$1=="xwasloc" { print $2 }' $paramfile ) # location of xwas executable
geneinterloc=$( awk '$1=="genescriptloc" { print $2 }' $paramfile ) # location of the gene_based interaction R script
genelist=$( awk '$1=="genelistname" { print $2 }' $paramfile ) # list of genes
genelist2=$( awk '$1=="genelistname2" { print $2 }' $paramfile ) # list of genes (optional)
buff=$( awk '$1=="buffer" { print $2 }' $paramfile ) # buffer around a gene (in base pairs)
covarfile=$( awk '$1=="covarfile" { print $2}' $paramfile) # covariate file (optional)
outputname=$( awk '$1=="output" { print $2 }' $paramfile ) # name of output file

echo -e "\tgene1\tgene2\tpmin\tpgate\tptail\tpprod" > ${outputname}

if [[ -z "${genelist2// }" ]]; then
  calc ${genelist} ${genelist}
else
  calc ${genelist} ${genelist2}
fi

echo "Wrote gene-gene interaction results to [ ${outputname} ]"

rm *.log
