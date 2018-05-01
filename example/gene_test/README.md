# Gene Based Association Test Example

To run the example, execute the shell script `run_example_gene_test.sh`.

To use XWAS gene-based tests, the following R packages must be installed,
`corpcor` and `mvtnorm`.

## Overview

This example demonstrates the use of ../../bin/gene_base_test.sh.
Gene-based testing builds upon the SNP-level analyses by using
the P -values obtained for each SNP.

## Requirements

The script for gene based testing requires:

1. A dataset in binary PLINK format
2. A file with a list of genes to test
3. A file with association statistics for SNPs
4. A file listing all the parameters necessary to run the script.

Example files with all of the above are provided for to in this example.

The output will be in `example_gene_test_result.txt`

## Reference

For more details please refer to the XWAS manual section 7, page 13.
