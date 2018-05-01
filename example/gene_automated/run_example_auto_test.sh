#!/bin/sh

if [ ! -f ./*_temp* ];
then
    ../../bin/xwas --noweb --bfile example2 --xwas --strat-sex --stouffers --out example2 && ../../bin/gene_based_test_automate.sh example_params_gene_test_auto.txt
else
    ../../bin/gene_based_test_automate.sh example_params_gene_test_auto.txt
fi

