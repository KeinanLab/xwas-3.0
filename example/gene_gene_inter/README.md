# Gene-based Gene-gene Interaction Testing example

To run the example, execute the shell script `run_example_GGG_test.sh`.

## Overview

This example demonstrates the use of ../../bin/gene_based_interaction.sh. This
test combines SNP-based interaction tests between all pairs of SNPs in two
genes to produce a gene-level test for interaction between the two.

Currently, the GGG tests can only deal with quantitative phenotypes. In order to
run the GGG tests, the following two packages must be installed in R: `corpcor`
and `mvtnorm`.

## Requirements

The script for gene based testing requires:

1. A dataset in binary PLINK format
2. A file with a list of genes to interaction for each pair of gene
3. A file listing all the parameters necessary to run the script.

Example files with all of the above are provided for to in this example.

The output will be in `genepair_result.txt`.

## Features

GGG test supports:
All X All, where All indicates all X chromosome genes.
SET1 X SET1, i.e. all pairs of genes in SET1
SET1 X SET2, i.e. all pairs of genes in SET1 and SET2
SET X ALL

Default behavior is to run on SET1 X SET1.

To run GGG on some cartesian product of all pairs of X chromosomes, gunzip
either the X chromosome gene list for either the hg19 or hg38 genome build in the
../../genes/ dir.

To run interaction tests on two differient sets of genes, add the optional
parameter `genelistname2` to the paramater file

```
filename    dummy_quant1
xwasloc	../../bin/
genescriptloc	../../bin/
genelistname	dummy_gene_test_list.txt
genelistname2	dummy_gene_test_list2.txt
buffer	0
output	genepair_result.txt
```

## Reference

For more details please refer to the XWAS manual section 8, page 14.
