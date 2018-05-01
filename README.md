XWAS
====

A Toolset for Chromosome X-Wide Association Studies

Built by the [Keinan Lab](http://keinanlab.cb.bscb.cornell.edu/) at Cornell University.

----
## Examples

Navigate to `$path/XWAS/example/sample_data` to find two sample datasets. Below are some commands you can use to run examples of XWAS tests:

1. Get allele frequencies for males and females separately
`../../bin/xwas --bfile dummy_case1 --xwas --freq-x`

2. Perform the sex-stratified test, using Fisher's method to combine p-values
`../../bin/xwas --bfile dummy_case1 --xwas --strat-sex --fishers`

3. Perform the sex-stratified test, using Stouffer's method to combine p-values
`../../bin/xwas --bfile dummy_case1 --xwas --strat-sex --stouffers`

4. Perform the sex-difference test
`../../bin/xwas --bfile dummy_quant1 --xwas --sex-diff`

5. Perform the variance-heterogeneity test (with covariates)
`../../bin/xwas --bfile dummy_quant1 --xwas --var-het --covar dummy_quant1.covar`

6. Perform linear regression using 0/2 coding for males
`../../bin/xwas --bfile dummy_quant1 --xchr-model 2 --linear`

For more information, [visit our website](http://keinanlab.cb.bscb.cornell.edu/content/xwas) to download the full XWAS manual.