#!/usr/bin/env python

import os
from sys import argv

# ARGUMENTS:
# 1. Old file root
# 2. New file root
old, new = argv[1], argv[2]
path = os.path.dirname(os.path.realpath(__file__))

# transform everything to binary datasets first
if not os.path.isfile(old + ".fam"):
    print("Converting " + old + " to binary dataset...")
    xwas_command = path + "/../bin/xwas --noweb --make-bed --file " + old + " --out " + old
    os.system(xwas_command)
if not os.path.isfile(new + ".fam"):
    print("Converting " + new + " to binary dataset...")
    xwas_command = path + "/../bin/xwas --noweb --make-bed --file " + new + " --out " + new
    os.system(xwas_command)

print("Comparing calls via XWAS...")
xwas_command = path + "/../bin/xwas --noweb --bfile " + old + " --bmerge " + new + ".bed " + new + ".bim " + new + ".fam --merge-mode 7 --out TEMP"
os.system(xwas_command)

print("Reading diff file...")
diff_file = "TEMP.diff"
snps = {}
individuals = {}
with open(diff_file) as diff:
    diff.next() # header
    for line in diff:
        line = line.split()
        snpid, iid = line[0], line[2]
        try:
            snps[snpid] += 1
        except KeyError:
            snps[snpid] = 1
        try:
            individuals[iid] += 1
        except KeyError:
            individuals[iid] = 1

old_snps = set()
old_individuals = set()
print("Reading " + old + " individuals...")
with open(old + ".fam") as old_fam:
    for line in old_fam:
        line = line.split()
        iid = line[1]
        old_individuals.add(iid)
print("Reading " + old + " snps...")
with open(old + ".bim") as old_bim:
    for line in old_bim:
        line = line.split()
        snpid = line[1]
        old_snps.add(snpid)

new_snps = set()
new_individuals = set()
print("Reading " + new + " individuals...")
with open(new + ".fam") as new_fam:
    for line in new_fam:
        line = line.split()
        iid = line[1]
        try:
            old_individuals.remove(iid)
        except KeyError:
            new_individuals.add(iid)
print("Reading " + new + " snps...")
with open(new + ".bim") as new_bim:
    for line in new_bim:
        line = line.split()
        snpid = line[1]
        try:
            old_snps.remove(snpid)
        except KeyError:
            new_snps.add(snpid)

incl_snps_diff = open("included_snps.diff","w")
print("Writing to included_snps.diff...")
print >> incl_snps_diff, '\t'.join(["SNP_ID","NUM_DIFF_CALLS"])
for snpid in snps.keys():
    print >> incl_snps_diff, '\t'.join([snpid,str(snps[snpid])])
incl_snps_diff.close()

exl_snps_diff = open("excluded_snps.diff","w")
print("Writing to excluded_snps.diff...")
print >> exl_snps_diff, '\t'.join(["SNP_ID",old,new])
for iid in old_snps:
    print >> exl_snps_diff, '\t'.join([iid,"1","0"])
for iid in new_snps:
    print >> exl_snps_diff, '\t'.join([iid,"0","1"])
exl_snps_diff.close()

incl_indiv_diff = open("included_individuals.diff","w")
print("Writing to included_individuals.diff...")
print >> incl_indiv_diff, '\t'.join(["IID","NUM_DIFF_CALLS"])
for iid in individuals.keys():
    print >> incl_indiv_diff, '\t'.join([iid,str(individuals[iid])])
incl_indiv_diff.close()

exl_indiv_diff = open("excluded_individuals.diff","w")
print("Writing to excluded_individuals.diff...")
print >> exl_indiv_diff, '\t'.join(["IID",old,new])
for iid in old_individuals:
    print >> exl_indiv_diff, '\t'.join([iid,"1","0"])
for iid in new_individuals:
    print >> exl_indiv_diff, '\t'.join([iid,"0","1"])
exl_indiv_diff.close()

os.system("rm TEMP*")
