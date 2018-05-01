# Automated Gene-based Association Tests

To run the example, execute, `run_example_auto_test.sh`.

Output will be in *_gene_results.txt

## Overview

To use XWAS gene-based tests, the following R packages must be installed,
`corpcor` and `mvtnorm`.

This example demonstrates the use of ../../bin/gene_based_test_automate.sh,
performing gene-based testing which builds upon previously performed
SNP-level analyses by using the P-values obtained for each SNP.

The script works by searching for the temporary files our XWAS software created
when you ran the SNP-level anlayses, and generates the gene based results for
you.

## Reference

For more details please refer to the XWAS manual section 7, page 12.
