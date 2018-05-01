#!/usr/bin/perl

use Getopt::Std;
use File::Basename;

$a=$ARGV[0];

### run smartpca
$parfile = "$a.par"; 
open(PARFILE,">$parfile") || die("OOPS couldn't open file $parfile for writing");
print PARFILE ("genotypename:	$a.ped\n");
print PARFILE ("snpname:	$a.map\n");
print PARFILE ("indivname:	$a.ped\n");
print PARFILE ("outputformat:	EIGENSTRAT\n");
print PARFILE ("genotypeoutname:	$a.eigenstratgeno\n");
print PARFILE ("snpoutname:	$a.snp\n");
print PARFILE ("indivoutname:	$a.ind\n");
print PARFILE ("familynames:	NO\n");

close(PARFILE);
$command = "./convertf";          # MUST put bin directory in path
$command .= " -p $parfile > $a.log";
system("$command");
