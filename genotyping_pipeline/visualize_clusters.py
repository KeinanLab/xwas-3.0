#!/usr/bin/env python2

import numpy as np
import csv
import argparse
import os
import matplotlib as mpl
import sys

try:
    os.environ["DISPLAY"]
except KeyError:
    print("WARNING: you are running without a $DISPLAY bash variable")
    print("You cannot view plots live, you must save them to files with -p")
    mpl.use('Agg')

import matplotlib.pyplot as plt

# Filter function to turn np sex array into bool array
# True values are female, False are male
def sex_filter(val):
    if val == 2: # female
        return True
    elif val == 1: # male
        return False
    else:
        return False

# Filter function to turn genotype array into color array
# Heterozygous are blue, ungenotyped are black
# Other homozygous are variety of colors
def genotype_filter(val):
    if val == 'AA':
        return 'cyan'
    elif val == 'CC':
        return 'red'
    elif val == 'GG':
        return 'magenta'
    elif val == 'TT':
        return 'green'
    elif val == '00': # ungenotyped
        return 'white'
    else: # heterozygous
        return 'blue'

# Joins a separated genotype array from tped
# [C, C, G, C] > [CC, GC]
def join_genotypes(sep_geno):
    join_geno = []
    i = 0
    while i != len(sep_geno):
        join_geno.append(sep_geno[i] + sep_geno[i+1])
        i += 2
    return join_geno

# Loads a list of SNP IDs from a text file into an array
def load_snps(snps):
    snp_set = set()
    with open(snps,'r') as f:
        for row in f:
            row = row.split()
            snp_set.add(row[0])
    return snp_set

# Load cel map used by Birdsuite into a dictionary
def load_cel_map(map):
    dictionary = {}
    with open(map,'r') as f:
        for row in f:
            row = row.split()
            cel = row[0]
            iid = row[1]
            dictionary[cel] = iid
    return dictionary

# Loads allele intensities from Birdsuite file for all SNPs in snps
# Returns A, B such that A[s][i] is the intensity of allele A in
# individual indiv_order[i] for SNP snp_order[s]
def load_intensities(allele_summary,snps):
    A = []
    B = []
    snp_order = []
    with open(allele_summary,'r') as f:
        row = f.next().split()
        while row[0].startswith("#"):
            row = f.next().split()
        header = row
        try:
            while len(snps) != 0:
                row = f.next().split()
                snpid = row[0][:-2]
                if snpid in snps:
                    A.append(row[1:])
                    B.append(f.next().split()[1:])
                    snps.remove(snpid)
                    snp_order.append(snpid)
        except StopIteration:
            sys.exit("ERROR: " + ', '.join(snps) + " do not exist in " + allele_summary)
    A = np.array(A).astype(float)
    B = np.array(B).astype(float)
    indiv_order = header[1:]
    return A,B,indiv_order,snp_order

# Load sex of individuals from .tfam file. Re-order according to
# allele_summary file order. Returns sex array and index mask
# to use on genotypes from .tped.
def load_tfam(tfile, celmap, indiv_order):
    # Make dictionary of what genotypes belong to what individuals
    dictionary = {} # k = IID, v = index in genotype array
    sex = []
    with open(tfile + ".tfam",'r') as f:
        index = 0
        for row in f:
            row = row.split()
            iid = row[1]
            dictionary[iid] = index
            index += 1
            sex.append(row[4])

    # Create an index re-ordering according to order of individuals in allele_summary
    if celmap:
        celdict = load_cel_map(celmap)

    index_map = []
    for iid in indiv_order:
        if iid[-4:] == ".CEL" or iid[-4:] == ".cel":
            iid = iid[:-4] # chop off '.CEL'
        try:
            if celmap:
                iid = celdict[iid]
            index_map.append(dictionary[iid])
        except KeyError:
            sys.exit("ERROR: " + iid + " does not exist in " + tfile + ".tfam")
    index_map = np.array(index_map)
    sex = np.array(sex)[index_map].astype(int)
    vfunc = np.vectorize(sex_filter)
    sex = vfunc(sex)
    return index_map,sex

# Load genotypes from .tped of all SNPs in snps. They should be presented
# in the same order they were encountered in .allele_summary
def load_tped(tfile,snps,snp_order):
    # Read in genotypes, order them like allele summary, and return
    genotypes = []
    index = 0
    with open(tfile + ".tped",'r') as f:
        try:
            while len(snps) != 0:
                row = f.next().split()
                snpid = row[1]
                if snpid in snps:
                    genotypes.append(join_genotypes(row[4:]))
                    snps.remove(snpid)
                    if snp_order[index] != snpid:
                        sys.exit("ERROR: " + tfile + " and allele_summary must have SNPs in the same order.")
                    index += 1
        except StopIteration:
            sys.exit("ERROR: " + ', '.join(snps) + " does not exist in " + tfile)
    vfunc = np.vectorize(genotype_filter)
    return vfunc(np.array(genotypes))

if __name__ == '__main__':

    BUFFER = 100

    parser = argparse.ArgumentParser(description='Arguments to visualize intensity clusters called by Birdsuite. See XWAS manual for further details. ')

    # Allele summary
    parser.add_argument('-a', '--allele-summary', action="store", dest="ALLELE_SUMMARY", type=str, required = True, help='.allele_summary file produced by Birdsuite')
    # TPED/TFAM file
    parser.add_argument('-t', '--tfile', action="store", dest="TFILE", type=str, required = True, help='a PLINK transposed fileset (.tped, .tfam) prefix')
    # SNP
    parser.add_argument('-s', '--snp', action="store", dest="SNPS", type=str, required = True, help='either a single SNP to plot to a text file with a list of SNPs')
    # CEL map
    parser.add_argument('-c', '--cel-map', action="store", dest="MAP", type=str, help='cel map used in Birdsuite')
    # Save PNG
    parser.add_argument('-p', '--png', action='store_true', dest='png', help = 'use this flag to save graph as a PNG')
    # Sex
    parser.add_argument('-f', '--female', action='store_true', dest='female', help = 'use this flag to view female clusters only')
    parser.add_argument('-m', '--male', action='store_true', dest='male', help = 'use this flag to view male clusters only')

    print("Parsing arguments...")
    args = parser.parse_args()

    if os.path.isfile(args.SNPS):
        print("Reading SNPs from file...")
        snps = load_snps(args.SNPS)
    else:
        print("Reading SNP name from command line...")
        snps = set([args.SNPS])

    print("Reading allele intensities...")
    A,B,indiv_order,snp_order = load_intensities(args.ALLELE_SUMMARY,snps.copy())
    print("Reading individuals...")
    index_map,sex = load_tfam(args.TFILE,args.MAP,indiv_order)
    print("Reading genotypes...")
    genotypes = load_tped(args.TFILE,snps,snp_order)

    if args.male:
        sex = np.invert(sex)

    for i in range(len(A)):
        print("Plotting " + snp_order[i] + "...")
        plt.xlim(0,max(A[i]) + BUFFER)
        plt.ylim(0,max(B[i]) + BUFFER)
        if args.female or args.male:
            plt.scatter(A[i][sex],B[i][sex],c=genotypes[i][index_map][sex],edgecolors='black')
            n = A[i][sex].size
        else:
            plt.scatter(A[i],B[i],c=genotypes[i][index_map],edgecolors='black')
            n = A[i].size
        if args.female:
            plt.title(snp_order[i] + ' females (n = ' + str(n) + ')')
            fname = snp_order[i] + '_female.png'
        elif args.male:
            plt.title(snp_order[i] + ' males (n = ' + str(n) + ')')
            fname = snp_order[i] + '_male.png'
        else:
            plt.title(snp_order[i] + ' (n = ' + str(n) + ')')
            fname = snp_order[i] + '.png'
        if args.png:
            plt.savefig(fname)
            plt.clf()
            plt.cla()
        else:
            plt.show()
