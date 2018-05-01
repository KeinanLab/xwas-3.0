#!/usr/bin/perl
use strict;
use warnings;

## This program is to check if SNPs in a bim file are concordant with reference SNPs from IMPUTE2;
## Date : 2016.01.21
## Programmer : yekaixiong

## Example useage: perl check_genome_build_and_strange_alignment.pl *.bim /path/to/reference_file.legend output_file_name.txt

## ^ is not a symbol in perl, rather I should use **
die "Usage:<input of bim file><input of reference file><output>" unless @ARGV == 3;

my $bim = $ARGV[0];
my $ref = $ARGV[1];
my @ref_files  = split /\//, $ref;
$ref_files[-1] =~ /chr([^_.]+)/;
my $chrname = $1;


open OA, ">$ARGV[2]" || die "OA:$!\n";
print OA "Chromosome $chrname is used for testing\n";

my %chr;
if($chrname eq 'X' || $chrname eq 'x'){
	$chr{23} = 1;
	$chr{'X'}= 1;
	$chr{'x'}= 1;
}else{
	$chr{$chrname} = 1;
}

open IA, "$bim"  || die "Can't open IA: $!\n";
my %hash;
my %SNPid;
my %counth;
my %countS;
my $removed = 0;
while (<IA>) {
	chomp;
	my @tmp = split;
   	next unless (defined $chr{$tmp[0]});
	if($tmp[4] eq '0' || $tmp[5] eq '0'){
		$removed++;
		next;
	}
	if(defined $hash{$tmp[3]}){
		print "Position $tmp[3] appear more than once\n";
		next;
	}else{
		my $alleles    = $tmp[4] . $tmp[5];
		$hash{$tmp[3]} = $alleles;
		$counth{$tmp[3]}= -9;
		my @inf = split /_/, $tmp[1];
		if($inf[0] =~ /^rs/){
			$SNPid{$inf[0]} = $alleles;
			$countS{$inf[0]}= -9;
		}
	}
}

close IA;
print OA "$removed SNPs removed because their alleles are 0\n";

open IB, "$ref" || die "$ref\n";
<IB>;
while(<IB>){
	chomp;
	my @tmp = split; 
	if(defined $tmp[4]){
		next unless ($tmp[4] eq "Biallelic_SNP");	
	}
	my @inf = split /:/, $tmp[0];
	my $bp1 = $tmp[2] . $tmp[3];
	my $bp2 = $tmp[3] . $tmp[2];
	if(defined $hash{$tmp[1]}){
		$counth{$tmp[1]} = 1;
		if ($hash{$tmp[1]} eq $bp1 || $hash{$tmp[1]} eq $bp2){
			$counth{$tmp[1]} = 2;
		}
	}

	next unless ($inf[0] =~ /^rs/);
	if(defined $SNPid{$inf[0]}){
		$countS{$inf[0]} = 1;
		if ($SNPid{$inf[0]} eq $bp1 || $SNPid{$inf[0]} eq $bp2){
			$countS{$inf[0]} = 2;
		}
	}

}
close IB;

print OA "\n";
print OA "Position-based concordance check:\n";
open OB, ">$ARGV[2].inconsistent" || die "OB:$!\n";
my $total_pos = 0;
my $found_pos = 0;
my $con_allele= 0;
foreach my $k (sort {$a<=>$b} keys %counth){
	$total_pos++;
	if($counth{$k} == 1){
		$found_pos++;
		print OB "$k\tPosition-based\n";
	}elsif($counth{$k} == 2){
		$found_pos++;
		$con_allele++;
	}
}

my $per1 = ($total_pos != 0)? (sprintf "%.2f", $found_pos/$total_pos * 100):('NA');
my $per2 = ($found_pos != 0)? (sprintf "%.2f", $con_allele/$found_pos*100) :('NA');

print OA "Total positions found in .bim: $total_pos\n";
print OA "Total positions found in .bim and in reference file: $found_pos ($per1%)\n";
print OA "Position-based allele concordance: $con_allele ($per2%)\n";
print OA "\n";

$total_pos = 0;
$found_pos = 0;
$con_allele= 0;
print OA "SNPid-based concordance check:\n";
foreach my $k (keys %countS){
        $total_pos++;
        if($countS{$k} == 1){
                $found_pos++;
                print OB "$k\tSNPid-based\n";
        }elsif($countS{$k} == 2){
                $found_pos++;
                $con_allele++;
        }
}

$per1 = ($total_pos != 0)? (sprintf "%.2f", $found_pos/$total_pos * 100):('NA');
$per2 = ($found_pos != 0)? (sprintf "%.2f", $con_allele/$found_pos*100):('NA');
print OA "Total IDs found in .bim: $total_pos\n";
print OA "Total IDs found in .bim and in reference file: $found_pos ($per1%)\n";
print OA "ID-based allele concordance: $con_allele ($per2%)\n";
