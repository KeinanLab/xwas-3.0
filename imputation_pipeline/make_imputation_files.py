#######################################################################################
#                 Generating Shell Scripts to do imputation                           #
#                    Feng Gao, Diana Chang and and Alexander Gorelick                 #
#                      Last update: September 22, 2015                                #
#                                                                                     #
# Usage: python PARFILE                                                               #
#                                                                                     #
# PARFILE has the following contents:                                                 #
#                                                                                     #
# FILE: name of pre-imputation PLINK files, WITHOUT EXTENSION!                        #
# OUTPUT: name of the output file                                                     #
# NJOBS: Number of parallel jobs                                                      #
# BUILD: Build of the data                                                            #
# FILELOC: location of the pre-imputation PLINK files                                 # 
# REFLOC: location of 1000G reference files                                           #
# TOOLSLOC: location of other imputation-related files                                #
# RESLOC1: intended location of the processed PLINK files (liftover, etc.)            #
# RESLOC2: intended location of the first round results from impute2 (before cat)     #
# FINALRESLOC: intended location of the final round results from impute2 (cat)        #
# MAFRULE: pop name (AFR.MAF, AMR.MAF, EAS.MAF, EUR.MAF, SAS.MAF, or ALL.MAF) and MAF #
# for filtering out the sites in the 1000G reference file                             #
#                                                                                     #
#######################################################################################

from os import system
from sys import argv
import os.path

fname, njob, bld, loc1, loc2, loc3, loc4, loc5, loc6, ffname, mafrule = [''] * 11

fin = open(argv[1], 'r')
lines = fin.readlines()
fin.close()
for line in lines:
	l = line.strip().split(':')
	if len(l) > 1:
		l1 = l[1].lstrip()
	if l[0] == 'FILE':
		fname = l1
	elif l[0] == 'OUTPUT':
		ffname = l1
	elif l[0] == 'NJOBS':
		njob = int(l1)
	elif l[0] == 'BUILD':
		bld = l1
# 	elif l[0] == 'FILTER':
# 		filt_rules = l1
	elif l[0] == 'FILELOC':
		loc1 = l1
	elif l[0] == 'REFLOC':
		loc2 = l1
	elif l[0] == 'TOOLSLOC':
		loc3 = l1
	elif l[0] == 'RESLOC1':
		loc4 = l1
	elif l[0] == 'RESLOC2':
		loc5 = l1
	elif l[0] == 'FINALRESLOC':
		loc6 = l1
	elif l[0] == 'MAFRULE':
		mafrule = l1

locFname = loc1 + fname
bed, bim, fam = locFname + '.bed', locFname + '.bim', locFname + '.fam'
geneticmap = loc2 + 'genetic_map_chrX_nonPAR_combined_b37.txt'
legend = loc2 + '1000GP_Phase3_chrX_NONPAR.legend'
#hapgz = loc2 + '1000GP_Phase3_chrX_NONPAR.hap.gz'
hap = loc2 + '1000GP_Phase3_chrX_NONPAR.hap'
if not (os.path.isfile(hap)):
	loc2 = '../'+loc2
	hap = loc2 + '1000GP_Phase3_chrX_NONPAR.hap'
	if not (os.path.isfile(hap)):
		print "At "+loc2+" nor "+loc2[3: ]+", cannot find 1000GP_Phase3_chrX_NONPAR.hap" #error msg
xwas, liftover, liftmefast, imputationsteps, gtool, impute2, catimputed = loc3 + 'xwas', loc3 + 'liftOver', loc3 + 'liftme_fast.sh', loc3 + 'imputationsteps.sh', loc3 + 'gtool', loc3 + 'impute2', loc3 + 'cat_imputed.sh'
hgchain = 'hg18ToHg19.over.chain' if bld == '18' else 'hg17ToHg19.over.chain'
hgchain = loc3 + hgchain
out1Name = fname + '_preimpute'  ## revised by kaixiong using "_" instead of "." per required by plink
pedgen, pedsample = out1Name + '.ped.gen', out1Name + '.ped.sample'
locpedgen, locpedsample = loc4 + pedgen, loc4 + pedsample 
out2Name = fname + '.impute2'

fout = open(fname + '_preimpute.sh', 'w')
fout.write('mkdir workdir\n')
fout.write('cp ' + bed + ' ' + bim + ' ' + fam + ' ' + xwas + ' ' + liftover + ' ' + liftmefast + ' ' + hgchain + ' ' + imputationsteps + ' ' + gtool + ' ./workdir\n')
fout.write('cd workdir\n')
fout.write('chmod 777 *.sh gtool liftOver xwas\n')
fout.write('./imputationsteps.sh ' + fname + ' ' + out1Name + ' ' + bld + '\n') 
fout.write('cd ../\n')	
fout.write('cp ./workdir/' + out1Name + '* ' + loc4 + '\n')
fout.write('rm -rf workdir\n')
fout.close()

system('chmod 777 ' + fname + '_preimpute.sh')


bg, ed = 2699520, 154931044 ## revised by Kaixiong as the coordinate of PAR1 & 2 in hg19
step = (ed - bg) / njob 
pts = [int(bg + step * i) for i in xrange(njob + 1)]
pts[njob] = ed ## revised by Lauren to round up last region to end of PAR
#print pts

for i in xrange(njob):
	fout = open(fname + '_impute2_run_' + `i` + '.sh', 'w')
	fout.write('mkdir workdir_' + `i` + '\n')
	fout.write('cp ' + impute2 + ' ' + loc4 + pedgen + ' ' + loc4 + pedsample + ' ./workdir_' + `i` + '\n')
	st = `pts[i]`
	if i == (njob - 1):
		fi = `pts[i + 1]`
	else:
		fi = `pts[i + 1] - 1`
	fout.write('cd workdir_' + `i` + '\n')	
	#fout.write('gzip -d 1000GP_Phase3_chrX_NONPAR.hap.gz\nchmod 777 impute2\n')
	fout.write('chmod 777 impute2\n')
	fout.write('./impute2 -chrX -m ../' + geneticmap + ' -h ../' + hap + ' -l ../' + legend + ' -g ' + pedgen + ' ' + '-sample_g ' + pedsample + ' -int ' + st + ' ' + fi + ' -Ne 20000 -o ' + out2Name + '.' + `i` + " -filt_rules_l 'TYPE!=Biallelic_SNP'" + " '" + mafrule + "'" + "\n")
	fout.write('cd ../\n')	
	fout.write('cp ./workdir_' + `i` + '/' + out2Name + '* ' + loc5 + '\n')
	fout.write('rm -rf workdir_' + `i` + '\n')	
	system('chmod 777 ' + fname + '_impute2_run_' + `i` + '.sh')

fout = open(fname + '_impute2_run_all.sh', 'w')
s = '\n'.join(['./' + fname + '_impute2_run_' + `i` + '.sh' for i in xrange(njob)])
fout.write(s)
fout.close()
system('chmod 777 ' + fname + '_impute2_run_all.sh')

fout = open(fname + '_impute2_cat.sh', 'w')
fout.write('mkdir catdir\n')
fout.write('cp ' + loc4 + out1Name + '.ped.sample ' + loc5 + '* ' + catimputed + ' ' + gtool + ' ' + xwas + ' ./catdir\n')
fout.write('cd catdir\n')
fout.write('chmod 777 *.sh gtool xwas\n')
fout.write('mv ' + out1Name + '.ped.sample ' + out2Name + '.ped.sample\n')
fout.write('./cat_imputed.sh ' + out2Name + ' ' + `njob` + ' 0.5 ' + out2Name + '.ped.sample ' + ffname + '\n')
fout.write('cd ../\n')
fout.write('cp ./catdir/*info ./catdir/' + ffname + '* ' + loc6 + '\n')

fout.write('rm -rf catdir\n')
fout.close()
system('chmod 777 ' + fname + '_impute2_cat.sh')	
