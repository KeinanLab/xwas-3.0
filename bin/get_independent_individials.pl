#!/usr/bin/perl
use strict;
use warnings;

## This program is to ;
## Date : 2015.
## Programmer : yekaixiong

## ^ is not a symbol in perl, rather I should use **
die "Usage:<input of base name for *.fam file>" unless @ARGV == 1 ;

if($ARGV[0] =~ /\.fam$/){
	die "Please input file basename only, without the extension of \"fam\"\n";
}

open OC, ">$ARGV[0].PC_Sib_clean.ind_to_remove" || die "IC:$!\n";

print "1: Removing individuals whose parents is/are also included.\n";
open IA, "$ARGV[0].fam"  || die "Can't open IA: $!\n";
open OA, ">$ARGV[0].PC_clean.removed" || die "OA:$!\n";
my %parents;
my %sample;
while (<IA>) {
	chomp;
	my @tmp = split;
   	$sample{$tmp[1]}  = 1;
	$parents{$tmp[2]} = 1 unless ($tmp[2] eq "0");
	$parents{$tmp[3]} = 1 unless ($tmp[3] eq "0");
}
close IA;
my %parents_as_sample;
my $num_p_as_s = 0;
foreach my $ind (keys %sample){
	if (defined $parents{$ind}){
		$parents_as_sample{$ind} = 1;
		$num_p_as_s++;
	}

}
print "$num_p_as_s parents are also included as sample. Their children should be removed.\n";
my $num_remove1 = 0;
open IA, "$ARGV[0].fam"  || die "Can't open IA: $!\n";
open OB, ">$ARGV[0].PC_clean.fam" || die "OA:$!\n";

my %times_as_parents;

while(<IA>){
	chomp;
	my @tmp = split; 
	if ( (defined $parents_as_sample{$tmp[2]}) || (defined $parents_as_sample{$tmp[3]})){
		my $tag;
		if(defined $parents_as_sample{$tmp[2]} && defined $parents_as_sample{$tmp[3]}){
			$tag = "Father and Mother Included as sample\n"
		}elsif(defined $parents_as_sample{$tmp[2]}){
			$tag = "Father Included as sample\n";
		}else{
			$tag = "Mother Included as sample\n";
		}
		print OA "$_\t$tag";	
		print OC "$tmp[0]\t$tmp[1]\n";
		$num_remove1++;
	}else{
		print OB "$_\n";
		$times_as_parents{$tmp[2]}++ unless ($tmp[2] eq "0");
		$times_as_parents{$tmp[3]}++ unless ($tmp[3] eq "0");
	}
}
print "$num_remove1 individuals are removed because their parent(s) is/are included as samples.\n";
close IA;
close OB;
close OA;
print "2: Removing individuals who have at least one siblings also included as samples (only first of the siblings in the file is retained).\n";
my $num_ind_with_siblings = 0;
my $num_ind_with_siblings_rved = 0;
my %retained;
open IB, "$ARGV[0].PC_clean.fam" || die "IB\n";
open OD, ">$ARGV[0].PC_Sib_clean.fam" || die "OD:$!\n";
open OE, ">$ARGV[0].Sib_clean.removed" || die "OA:$!\n";

$times_as_parents{0} = 1; ## This setting is to skip the analysis for samples with missing parental information

while(<IB>){
	chomp;
	my @tmp = split; 
	if ( ($times_as_parents{$tmp[2]} != 1) || ($times_as_parents{$tmp[3]} != 1)){
		$num_ind_with_siblings++;		
		if( (defined $retained{$tmp[2]}) || (defined $retained{$tmp[3]})){ ## this sample has one sibling included, so it has to be removed
			$num_ind_with_siblings_rved++;
			my $tag1 = (defined $retained{$tmp[2]})? $retained{$tmp[2]}:"NA";
			my $tag2 = (defined $retained{$tmp[3]})? $retained{$tmp[3]}:"NA";
			print OE "$_\t$tag1\t$tag2\n";
			print OC "$tmp[0]\t$tmp[1]\n";
		}else{ ## this sample has siblings, but it appears first in the file and will be included
			print OD "$_\n";
			$retained{$tmp[2]} = $tmp[1];
			$retained{$tmp[3]} = $tmp[1];
		}		

	}else{
		print OD "$_\n";
	}

}
print "$num_ind_with_siblings individuals have siblings in the sample.\n";
print "$num_ind_with_siblings_rved individuals removed as a result.\n";
close IB;
close OC;
close OD;
close OE;
