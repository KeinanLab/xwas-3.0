#!/usr/bin/perl
use strict;
use warnings;

## This program is to ;
## Date : 2015.
## Programmer : yekaixiong

## ^ is not a symbol in perl, rather I should use **
die "Usage:<input><output>" unless @ARGV == 2;


open IA, "$ARGV[0]"  || die "Can't open IA: $!\n";
open OA, ">$ARGV[1]" || die "OA:$!\n";
my %hash;
while (<IA>) {
	chomp;
	my @tmp = split;
   	$tmp[2] = 0;
	$tmp[3] = 0;
	my $line= join "\t", @tmp;
	print OA "$line\n";

}

close IA;
close OA;
