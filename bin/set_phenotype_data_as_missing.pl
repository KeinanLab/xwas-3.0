#!/usr/bin/perl
use strict;
use warnings;

## This program is to set phenotypic values as missing in *fam file;
## Date : 2015.01.01
## Programmer : yekaixiong

## ^ is not a symbol in perl, rather I should use **
die "Usage:<input of fam file to revise>" unless @ARGV == 1 ;


open IA, "$ARGV[0]"  || die "Can't open IA: $!\n";
open OA, ">new.fam" || die "OA:$!\n";
my $num = 0;
while (<IA>) {
	chomp;
	my @tmp = split;
	$num++ if ($tmp[-1] != -9);
	$tmp[-1]= -9;
	my $line= join "\t", @tmp;
	print OA "$line\n";
}
close IA;
close OA;
print "$num individuals have non-missing phenotype and were set as -9\n";
`mv new.fam $ARGV[0]`;
