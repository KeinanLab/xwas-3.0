# Quality Control step example

To run the example, execute, `run_example_qc.sh`.

In addition the parameter file, a few select flags can be directly passed to the
QC script by editing `run_example_qc.sh`.
* `-l, --save-logs`           saves logs from QC procedure to
         ./[FILENAME]_QC_logs
* `-a, --save-all`            saves all intermediate files from QC
         procedure, including logs, to ./[FILENAME]_QC_intermediate_files
* `-v, --verbose`             unsupresses plink output
* `-g, --debug`               saves logs, intermediate files, and keep plink output

( i.e. to run the QC step in debug mode `../../bin/run_QC.sh --debug
temp_params.txt` )

## Overview

In this example we demonstrate the usage of our pre-quality control and quality
control pipeline that applies standard autosomal GWAS quality control steps as
implemented in PLINK (Purcell et al. 2007) and SMARTPCA (Price et al. 2006), as
well as procedures that are specific to X. It is recommended (not required) to
perform imputation on the dataset using IMPUTE2

## Reference

For more details please refer to the XWAS manual section 3, page 4.
