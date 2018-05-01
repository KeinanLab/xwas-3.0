#!/usr/bin/perl
use strict;
use warnings;

## This program is to ;
## Date : 2015.12.31
## Programmer : yekaixiong

## ^ is not a symbol in perl, rather I should use **
die "Usage:<input><output of ind to remove>" unless @ARGV == 2;


open IA, "$ARGV[0]"  || die "Can't open IA: $!\n";
open OA, ">$ARGV[1]" || die "OA:$!\n";
<IA>;
my %hash;
my %ind;
my $num_pairs = 0;
my %family;
while (<IA>) {
	chomp;
	my @tmp = split;
   	$hash{$tmp[1]}->{$tmp[3]} = 1;
	$ind{$tmp[1]}++;
	$ind{$tmp[3]}++;
	$family{$tmp[1]} = $tmp[0];
	$family{$tmp[3]} = $tmp[2];
	$num_pairs++;
}
close IA;
#print "Total number of related pairs: $num_pairs\n";
my %removed;
my @sorted_ind = sort {$ind{$b} <=> $ind{$a}} keys %ind;
unless(defined $sorted_ind[0]){
	exit;
}
my $rmvind     = $sorted_ind[0];
my $rmpair     = $ind{$rmvind};

while($num_pairs != 0){
	$removed{$rmvind} = $rmpair;
	#print "Removed individual $rmvind, involving in $rmpair pairs\n";
	$num_pairs = 0;
	my %new_ind;
	foreach my $k1 (keys %hash){
		delete $hash{$k1} if (defined $removed{$k1});
		foreach my $k2 (keys %{$hash{$k1}}){
			if (defined $removed{$k2}){
				delete $hash{$k1}->{$k2};
			}else{
				$num_pairs++;
				$new_ind{$k1}++;
				$new_ind{$k2}++;
			}
		}
	}
	#print "Current number of related pairs: $num_pairs\n";
	if($num_pairs > 0){
		my @new_sorted_ind = sort {$new_ind{$b} <=> $new_ind{$a}} keys %new_ind;
		$rmvind = $new_sorted_ind[0];
		$rmpair = $new_ind{$new_sorted_ind[0]};
	}
}

foreach my $k (keys %removed){
	print OA "$family{$k}\t$k\trelated\n";
}
close OA;
