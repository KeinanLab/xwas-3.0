#!/usr/bin/perl
use strict;
use warnings;

## This program is to compare SNPs in male and female, exporting a list of SNPs not shared between the two. These SNPs should be removed before merging the male and female dataset;
## Date : 2016.1.18
## Programmer : yekaixiong

## ^ is not a symbol in perl, rather I should use **
die "Usage:<input of male bim file><input of female bim file><output of un-shared SNPs>" unless @ARGV == 3;


open IA, "$ARGV[0]"  || die "Can't open IA: $!\n";
my %hash;
while (<IA>) {
	chomp;
	my @tmp = split;
	$hash{$tmp[1]}++;
}
close IA;

open IB, "$ARGV[1]" || die "IB:$!\n";
while(<IB>){
	chomp;
	my @tmp = split; 
	$hash{$tmp[1]}++
}
close IB;
open OA, ">$ARGV[2]" || die "OA:$!\n";
foreach my $k (keys %hash){
	if($hash{$k} == 1){
		print OA "$k\n";
	}else{
		die "$k\t$hash{$k}\n" unless ($hash{$k} == 2);
	}
}
close OA;
