
///////////////////////////////////////////////////////////////////////////
//       Adapted from plink source codes by Purcell et al. 2007          //
//       By Feng Gao and Diana Chang                                     //
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <iterator>
#include <algorithm>
#include "options.h"
#include "helper.h"
#include "stats.h"
#include "crandom.h"
#include "linear.h"
#include "logistic.h"
#include "fisher.h"
#include "plink.h"
#include "whap.h"
#include "phase.h"
#include "xoptions.h"
#include "dcdflib.h"
#include <cfloat>

#define MISSING(i,l) ( P.SNP[l]->one[i] && ( ! P.SNP[l]->two[i] ) )

void xfilterSNPs(Plink& P)
{ 
	//////////////////////////////////////////////////////
	// This functions applies the following filters and
	// functions:
	// Counts of number of founders, nonfounders
	// Per-individual genotyping rate
	// Read in, or calculate, allele frequencies (and save these)
	// Optionally write allele frequencies, then close
	// Exclude SNPs with too many missing genotypes
	// Identify/correct heterozygote haploid
	// Identify SNPs with no founder genotypes
	// Calculate/report/filter on HWE tests
	// Calculate/report genotyping rate per SNP/per individual
	// Filter on MAF
	// Remove filtered-out SNPs


	bool original_SNP_major = par::SNP_major;

	if ( ! par::SNP_major ) P.Ind2SNP();

	// Which SNPs to delete
	vector<bool> del(P.locus.size(),false);

	// Which individuals to delete
 	vector<bool> indel(P.sample.size(),false);

	//////////////////////////////////////////
	// Display number of founders/nonfounders
	P.cnt_f=0;
	vector<Individual*>::iterator person = P.sample.begin();
	while ( person != P.sample.end() )
	{
		if ( (*person)->founder ) P.cnt_f++;
		person++;
	}
	P.printLOG(int2str(P.cnt_f)+" founders and "+int2str(P.n-P.cnt_f)+" non-founders found\n");
	if (P.cnt_f<P.n) par::has_nonfounders = true;
 
	////////////////////////////////////////////
	// If we have obligatory missing genotypes:
	// ensure they really are missng
	if ( par::oblig_missing ){
		set<int2>::iterator p = P.oblig_missing.begin();
		while ( p != P.oblig_missing.end() ){
			int l = p->p1;
			int k = p->p2;
			for (int i = 0; i < P.sample.size(); i++){
				Individual * person = P.sample[i];
				if ( person->sol == k ){
					P.SNP[l]->one[i] = true;
					P.SNP[l]->two[i] = false;
				}
			}
			++p;
		}
	}

	/////////////////////////////////////////////////
	// Remove individuals with too many missing calls

	double total_genotyping = 0;

	if ( par::MAX_IND_MISSING < 1 ){
		int n_removed = 0;
		int n_orig = P.n;
		// Consider each individual
		if ( ! par::oblig_missing ){
			for (int i = 0;i < P.sample.size(); i ++){
				bool female = !(P.sample[i]->sex);
				// Sum missingness over all SNPs
				int m=0;       // Missing SNPs
				int nsnps=0;   // All non-obligatory missing SNPs
				for (int l = 0; l < P.locus.size(); l ++){
					// Skip female Y chromosomes
					if ( female && par::chr_Y[P.locus[l]->chr] ) continue;
		  			++nsnps;
					if ( MISSING(i,l) ) m++;
				}
				// Too much missingness?
				if ( (double)m/(double)nsnps > par::MAX_IND_MISSING ){
					indel[i] = true;
					n_removed++;
				}
			} // next individual
		}
		else{ // ... allow oblig missing values
			for (int i=0;i<P.sample.size();i++)
			{
				bool female = ! (P.sample[i]->sex);
				// Sum missingness over all SNPs
				int m=0;       // Missing SNPs
				int nsnps=0;   // All non-obligatory missing SNPs
				for (int l=0; l<P.locus.size();l++){
		  			// Skip female Y chromosomes
					if (female && par::chr_Y[P.locus[l]->chr]) continue;
					if ( ! (P.obligMissing(i,l)) ){
						if (MISSING(i,l)) ++m;
						++nsnps;
					}
				}
				// Too much missingness?
				if ( (double)m/(double)nsnps > par::MAX_IND_MISSING ){
					indel[i] = true;
					n_removed++;
				}
			} // next individual
		} // end if oblig-missing section
		
		////////////////////////////////////////
		// Save list of any removed individuals
		
		if (n_removed>0){
			string f = par::output_file_name + ".irem";
			P.printLOG("Writing list of removed individuals to [ " + f + " ]\n");
			ofstream REM;
			REM.open(f.c_str(), ifstream::out);
			for (int i=0;i<P.sample.size();i++)
				if (indel[i]) REM << P.sample[i]->fid << "\t" << P.sample[i]->iid << "\n";
			REM.close();
			
			// And now remove these individuals, so that
			// SNP-based statistics are calculated with
			// these samples already excluded
			n_removed = P.deleteIndividuals(indel);
		}

		P.printLOG(int2str(n_removed)+" of "+int2str(n_orig));
		P.printLOG(" individuals removed for low genotyping ( MIND > ");
		P.printLOG(dbl2str(par::MAX_IND_MISSING)+" )\n");
	} // end of remove people conditional


	/////////////////////////////////
	// Calculate or read from file?

	if (par::af_read){
		checkFileExists(par::af_file);
		P.printLOG( "Reading allele frequencies from [ " + par::af_file + " ] \n");

		// Make hash of original SNP names
		map<string,int> mlocus;
		map<string,int>::iterator ilocus;

		vector<Locus*>::iterator loc = P.locus.begin();
		int l=0;
		while ( loc != P.locus.end() ){
			mlocus.insert(make_pair( (*loc)->name,l));
			loc++;
			l++;
		}

		// Read allele frequencies
		ifstream FRQ;
		FRQ.open(par::af_file.c_str());
		FRQ.clear();

		string dum1, dum2, dum3, dum4, dum5, dum6;
		string snpname;
		double freq;
		int nm;
		loc = P.locus.begin();
		while ( loc != P.locus.end() ){
			(*loc)->freq = -1;
			(*loc)->nm = 0;
			loc++;
		}
		// Skip header line
		FRQ >> dum1 >> dum2 >> dum3 >> dum4 >> dum5 >> dum6;
		while(!FRQ.eof()){
			vector<string> tokens = tokenizeLine(FRQ);
			if (tokens.size() == 0) continue;
			else if (tokens.size() != 6){
				string sline="";
				for (int i=0; i<tokens.size(); i++)
					sline += tokens[i] + " ";
					error("Problem with allele frequency line: 6 fields required:\n" + sline + "\n");
			}
			ilocus = mlocus.find(tokens[1]);
			if (ilocus != mlocus.end()){
				string rareAllele = tokens[2];
				string commonAllele = tokens[3];
				Locus * loc = P.locus[ilocus->second]; //not very sure
				if( ! from_string<double>( loc->freq, tokens[4],std::dec)){
					loc->freq = 0;
					loc->nm = 0;
				}
				else if( ! from_string<int>(loc->nm,tokens[5],std::dec)){
					loc->freq = 0;
					loc->nm = 0;
				}
				// But was that pointing to the correct allele?
				if ( rareAllele == loc->allele2 && rareAllele != par::missing_genotype && loc->allele2 != par::missing_genotype )
					loc->freq = 1 - loc->freq;
				else if ( commonAllele == loc->allele1 && commonAllele != par::missing_genotype && loc->allele1 != par::missing_genotype )
					loc->freq = 1 - loc->freq;
			}
		}
		FRQ.clear();
		FRQ.close();
	}

	/////////////////////////////////
	// Calculate allele frequencies

	vector<string> hetlist(0);
	vector<bool>::iterator d = del.begin();
	vector<Locus*>::iterator loc = P.locus.begin();
	vector<CSNP*>::iterator s = P.SNP.begin();
	int l = 0; // Main locus counter
	int exc_maf = 0;
	int exc_miss = 0;

	vector<Locus*> no_founders_found_list;
	vector< double > freq_f(P.locus.size(), 0), freq_m(P.locus.size(), 0); //added by Feng Gao, Diana Chang
	vector< double > miss_f(P.locus.size(), 0), miss_m(P.locus.size(), 0); //counter for missingness
	vector< int > nm_f(P.locus.size(), 0), nm_m(P.locus.size(), 0); //added by Feng Gao, Diana Chang, for allele frequencies of chrX for males and females separately.
	vector< string > allele_1(P.locus.size()), allele_2(P.locus.size()); //added by Feng Gao & Diana Chang. Record the two alleles when the freq calculation is done to avoid the swap of alleles.

	vector< string >::iterator a1_it = allele_1.begin(); //added by Feng Gao, Diana Chang
	vector< string >::iterator a2_it = allele_2.begin(); //added by Feng Gao, Diana Chang
	vector< double >::iterator freq_f_it = freq_f.begin(); //added by Feng Gao, Diana Chang
	vector< double >::iterator freq_m_it = freq_m.begin(); //added by Feng Gao, Diana Chang
	vector< double >::iterator miss_f_it = miss_f.begin(); //added by Arjun Biddanda (missingness in females per locus)
	vector< double >::iterator miss_m_it = miss_m.begin(); // added by Arjun Biddanda (missingness in males per locus)

	vector< int >::iterator nm_f_it = nm_f.begin(); //added by Feng Gao, Diana Chang
	vector< int >::iterator nm_m_it = nm_m.begin(); //added by Feng Gao, Diana Chang
	int exc_sexdiff = 0; //added by Feng Gao, Diana Chang
	int exc_missdiff = 0; //added by Arjun Biddanda
	while ( loc != P.locus.end() ){//go over each loci
		if (!par::af_read){
			(*loc)->freq = 0;
			// count 1 per allele, for frequency
			(*loc)->nm = 0;
			*freq_f_it = 0; //added by Feng Gao, Diana Chang
			*freq_m_it = 0; //added by Feng Gao, Diana Chang
			*miss_f_it = 0; //added by Arjun Biddanda
			*miss_m_it = 0; // added by Arjun Biddanda
			*nm_f_it = 0; //added by Feng Gao, Diana Chang
			*nm_m_it = 0; //added by Feng Gao, Diana Chang
		}
		*a1_it = (*loc)->allele1; //added by Feng Gao & Diana Chang
		*a2_it = (*loc)->allele2; //added by Feng Gao & Diana Chang
		// count 1 per genotype, for missingness
		int geno_nm = 0;

		// count 1 per non-obligatory missing genotype
		// (or set to N individuals)
		int geno_real = 0;
		bool X = false;
		bool haploid = false;

		// Determine type of SNP
		if (par::chr_sex[(*loc)->chr]) X=true; //Is this loc a X chr?
		else if (par::chr_haploid[(*loc)->chr]) haploid=true;

		///////////////////////////////
		// Iterate over each individual

		vector<bool>::iterator i1 = (*s)->one.begin();
		vector<bool>::iterator i2 = (*s)->two.begin();
		vector<Individual*>::iterator person = P.sample.begin();

		int i = 0;
		while ( person != P.sample.end() ){	
			bool s1 = *i1;
			bool s2 = *i2;
			// Check female Y genotypes
			if ( par::chr_Y[(*loc)->chr] && ! (*person)->sex ){
				// Set to missing, unless in a RECODE mode
				if ( ! par::preserve_all_genotypes ){
					s1 = *i1 = true;
					s2 = *i2 = false;
				}
				// But in any case, do not include this marker in
				// any genotype counts: skip to next person
				++person;
				++i;
				++i1;
				++i2;
				continue;
			}

			// For haploid heterozygosity check, also consider all individuals
			if ( haploid || ( X && (*person)->sex ) ){
				if ( (!s1) && s2 ){
					hetlist.push_back( (*person)->fid + "\t" + (*person)->iid + "\t" + (*loc)->name );
					// Set to missing, unless in a RECODE mode
					if ( ! par::preserve_all_genotypes ){
						s1 = *i1 = true;
						s2 = *i2 = false;
					}
				}	
			}
		
			// For missing genotypes (Arjun Biddanda)
			if ( !(s1 && (!s2)) ){
			 	geno_nm++;
			 	// Only calculate missingness in controls
			 	if (!(*person)->aff){
			 		if ((*person)->sex){
			 			*miss_m_it += 1;
			 		}
			 		else{
			 			*miss_f_it += 2;
			 		}
			 	}
			}
			// But is this a real genotype in any case?
			if ( par::oblig_missing ){
				if ( ! (P.obligMissing(i,l)) ) ++geno_real;
			}
			else ++geno_real;
			 



			// Do not recount alleles if we have read in allele frequencies
			if (!par::af_read){
				// For allele frequencies
				// only consider founders
				if ( par::summ_nonfounders || (*person)->founder ){
					// haploid or male(1) on X
					if ( haploid || ( X && (*person)->sex ) ){//chrX for males should be treated differently
						//////////////////
						// Haploid counts
						// "1" allele count

						if ( (!s1) && (!s2) ){   //  FF = hom(11), false = allele1 according to input.cpp
							(*loc)->freq++;
							(*loc)->nm++;
							if (!(*person)->aff) { //added by Feng Gao & Diana Chang 
								*freq_m_it += 1;
								*nm_m_it += 1;
							}
						}
						else if ( s1 && s2 ){   //  TT = hom(22)
							(*loc)->nm++;
							if (!(*person)->aff) { //added by Feng Gao & Diana Chang
								*nm_m_it += 1;
							}
						}
					}
					else{
						//////////////////
						// Autosomal count 
						// "1" allele count

						if (!s1){
							if (!s2){ //   00 = hom(11)
								(*loc)->freq += 2;
								(*loc)->nm += 2;
								if ((*person)->sex) { // added by Feng Gao & Diana Chang, unaffected male
									if (!(*person)->aff) {
										*freq_m_it += 2;
										*nm_m_it += 2;
									}
								} 
								else { // added by Feng Gao & Diana Chang, female
									if (!(*person)->aff) { 
										*freq_f_it += 2;
										*nm_f_it += 2;
									}
								}
							}
							else{ //   01 = het(12)
								(*loc)->freq += 1;
								(*loc)->nm += 2;
								if ((*person)->sex){ // added by Feng Gao & Diana Chang
									if(!(*person)->aff) { 
										*freq_m_it += 1;
										*nm_m_it += 2;
									}
								} 
								else { // added by Feng Gao & Diana Chang
									if(!(*person)->aff) { // unaffected
										*freq_f_it += 1;
										*nm_f_it += 2;
									}
								}
							}
						}
						else if ( s2 ){ // 11 = hom(22)
							(*loc)->nm+=2;
							if ((*person)->sex && !(*person)->aff){ // added by Feng Gao & Diana Chang
								*nm_m_it += 2;
							}
							else if (!(*person)->sex && !(*person)->aff){ // added by Feng Gao & Diana Chang
								*nm_f_it += 2;
							}
						}
					}
				}
			}

			// Next individual
			++person;
			++i;
			++i1;
			++i2;
		}


		////////////////////////////////
		// Calculate allele frequencies

		if (!par::af_read){
			if ((*loc)->nm>0) (*loc)->freq /= (double)(*loc)->nm;
			else{
				(*loc)->freq = 1;
				// If we aren't getting rid of it anyway
				if ( (double)geno_nm/(double)geno_real >= (1-par::MAX_GENO_MISSING)) no_founders_found_list.push_back(*loc);
			}
		}
		
		//////////////////////////////////////////
		// added by Feng Gao & Diana Chang. Calculating sex stratified X AF
		if (xpar::af_sex) {
			// Females
 			if (*nm_f_it > 0) *freq_f_it /= (double)*nm_f_it;
			else *freq_f_it = 1;
			// Males
			if (*nm_m_it > 0) *freq_m_it /= (double)*nm_m_it;
			else *freq_m_it = 1;
		}

		//added by Arjun Biddanda. Calculating differential missingness on X b/w M/F
		if (xpar::am_sex_test){

			//Females
			if (*miss_f_it > 0) *miss_f_it /= (double)*nm_f_it;
			else *miss_f_it = 1;
			//Males
			if (*miss_m_it > 0) *miss_m_it /= (double)*nm_m_it;
			else *miss_m_it = 1;

		}


		//////////////////////////////////////////
		// Record total proportion of missingness

		double snp_genotyping = P.n>0 ? (double)geno_nm/(double)geno_real : 0;
		if (snp_genotyping < 1.0){
			// cout<"SOME MISSING!\n";
		}
		total_genotyping += snp_genotyping;

		/////////////////////////////////////////////////
		// Exclude if SNP has too many missing genotypes
		if ( snp_genotyping < (1-par::MAX_GENO_MISSING) ){
			*d = true;
			exc_miss++;
		}

		////////////////////////////////////////////////
		// Make allele1 always the least common allele

		if ( par::make_minor_allele && (!par::af_count) && (*loc)->freq > 0.5 ){

			// then we need to swap alleles
			(*loc)->freq = 1 - (*loc)->freq;			
			string tmp = (*loc)->allele2;
			(*loc)->allele2 = (*loc)->allele1;
			(*loc)->allele1 = tmp;

			vector<bool>::iterator i1 = (*s)->one.begin();
			vector<bool>::iterator i2 = (*s)->two.begin();

			while ( i1 != (*s)->one.end() ){
				if ( (*i1) == (*i2) ){
					*i1 = ! (*i1);
					*i2 = ! (*i2);
				}
				i1++;
				i2++;
			}
		}

		// Next SNP
		++d;
		++loc;
		++l;
		++s;
		++a1_it;
		++a2_it;
		++freq_f_it;
		++freq_m_it;
		++miss_f_it;
		++miss_m_it;
		++nm_f_it;
		++nm_m_it;
	}

	/////////////////////////////////////////////////
	// Save list of any heterozygous haploid alleles

	if (hetlist.size()>0){
		P.printLOG(int2str( hetlist.size()) + " heterozygous haploid genotypes; set to missing\n");
		string f = par::output_file_name + ".hh";
		P.printLOG("Writing list of heterozygous haploid genotypes to [ " + f + " ]\n");
		ofstream REM;
		REM.open(f.c_str(), ifstream::out);
		for (int i=0; i<hetlist.size(); i++) REM << hetlist[i] << "\n";
		REM.close();
	}
	hetlist.clear();

	/////////////////////////////////////////////////
	// Save list of SNPs with no founders observed

	if (no_founders_found_list.size()>0){
		P.printLOG(int2str( no_founders_found_list.size()) + " SNPs with no founder genotypes observed\n");
		P.printLOG("Warning, MAF set to 0 for these SNPs (see --nonfounders)\n");
		string f = par::output_file_name + ".nof";
		P.printLOG( "Writing list of these SNPs to [ " + f + " ]\n");
		ofstream NOF;
		NOF.open(f.c_str(), ifstream::out);
		for (int i=0; i<no_founders_found_list.size(); i++) NOF << no_founders_found_list[i]->name << "\n";
		NOF.close();
	}
	no_founders_found_list.clear();

	//////////////////////////
	// Write allele freq file

	if (par::af_write){
		if (par::include_cluster_from_file) P.calcStratifiedAlleleFreqs();
		else{
			ofstream FRQ;
			string f = par::output_file_name + ".frq";
			if (par::af_count) f += ".count";
			if (par::summ_nonfounders) P.printLOG("Writing allele frequencies (all individuals) to [ " + f + " ] \n");
			else P.printLOG("Writing allele frequencies (founders-only) to [ " + f + " ] \n");
			if (par::af_count) P.printLOG("Display counts rather than frequencies\n");
			FRQ.open(f.c_str(), ifstream::out);
			FRQ.precision(4);
			FRQ << setw(4) << "CHR" << " "
			<< setw(par::pp_maxsnp) << "SNP" << " "
			<< setw(4) << "A1" << " "
			<< setw(4) << "A2" << " ";
			if (par::af_count)
				FRQ << setw(6) << "C1" << " "
				<< setw(6) << "C2" << " "
				<< setw(6) << "G0" << "\n";
			else
				FRQ << setw(12) << "MAF" << " "
				<< setw(8) << "NCHROBS"
				<< "\n";

			vector<Locus*>::iterator loc = P.locus.begin();
			while (loc != P.locus.end() ){
				string a1 = (*loc)->allele1;
				string a2 = (*loc)->allele2;
				if (a1=="") a1="0";
				if (a2=="") a2="0";
				FRQ << setw(4)  << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp) << (*loc)->name  << " "
				<< setw(4)  << a1  << " "
				<< setw(4)  << a2  << " ";

				if (par::af_count){
					FRQ << setw(6) << int( (*loc)->freq ) << " "
					<< setw(6) << int( (*loc)->bp   ) << " "
					<< setw(6) << int( (*loc)->pos  ) << "\n";
				}
				else{
					if ( (*loc)->nm > 0 )	FRQ << setw(12) << (*loc)->freq << " ";
					else	FRQ << setw(12) << "NA" << " ";
					FRQ << setw(8) << (*loc)->nm << "\n";
				}
				loc++;
			}
			FRQ.close();
		}

		// Close after we've done alle freqs,
		shutdown();
	}

	//////////////////////////
	// added by Feng Gao & Diana Chang. Write allele freq file for X alleles

	if (xpar::af_sex && !par::af_read){
		ofstream FRQ;
		string f = par::output_file_name + ".xfrq";
		if (par::summ_nonfounders)	P.printLOG("Writing X sex stratified allele frequencies (all individuals) to [ " + f + " ] \n");
		else	P.printLOG("Writing X sex stratified allele frequencies (founders-only) to [ " + f + " ] \n");
		FRQ.open(f.c_str(), ifstream::out);
		FRQ.precision(4);
		FRQ << setw(4) << "CHR" << " "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(4) << "A1_F" << " "
		<< setw(4) << "A2_F" << " "
		<< setw(12) << "MAF_F" << " "
		<< setw(8) << "NCHROBS_F" << " "
		<< setw(4) << "A1_M" << " "
		<< setw(4) << "A2_M" << " "
		<< setw(12) << "MAF_M" << " "
		<< setw(8) << "NCHROBS_M" << "\n";
		vector<Locus*>::iterator loc = P.locus.begin();
		vector< string >::iterator a1_it = allele_1.begin();
		vector< string >::iterator a2_it = allele_2.begin();
		vector< double >::iterator freq_f_it = freq_f.begin(); 
		vector< double >::iterator freq_m_it = freq_m.begin(); 
		vector< int >::iterator nm_f_it = nm_f.begin(); 
		vector< int >::iterator nm_m_it = nm_m.begin(); 
		while (loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){ //only if this locus is on chrX
				string a1 = *a1_it;
				string a2 = *a2_it;
				if (a1=="") a1="0";
				if (a2=="") a2="0";
				string a1_f = (*freq_f_it <= (1. - *freq_f_it)) ? a1 : a2;
				string a2_f = (*freq_f_it <= (1. - *freq_f_it)) ? a2 : a1;
				string a1_m = (*freq_m_it <= (1. - *freq_m_it)) ? a1 : a2;
				string a2_m = (*freq_m_it <= (1. - *freq_m_it)) ? a2 : a1;
				FRQ << setw(4)  << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp) << (*loc)->name  << " "
				<< setw(4)  << a1_f  << " "
				<< setw(4)  << a2_f << " ";
				if ( *nm_f_it > 0 )	FRQ << setw(12) << min(*freq_f_it, 1. - *freq_f_it) << " ";
				else	FRQ << setw(12) << "NA" << " ";
				FRQ << setw(8) << *nm_f_it << " "
				<< setw(4)  << a1_m  << " "
				<< setw(4)  << a2_m << " ";
				if ( *nm_m_it > 0 )	FRQ << setw(12) << min(*freq_m_it, 1. - *freq_m_it) << " ";
				else	FRQ << setw(12) << "NA" << " ";
				FRQ << setw(8) << *nm_m_it << "\n";
			}
			loc ++;
			a1_it ++;
			a2_it ++;
			nm_f_it ++;
			nm_m_it ++;
			freq_f_it ++;
			freq_m_it ++;
		}
		FRQ.close();
		// Close after we've done alle freqs,
		shutdown();
	}

	//////////////////////////
	// added by Feng Gao & Diana Chang.
	// Check for significant differences between male / female allele counts
	// using fisher's exact test, then write results to a file.
	// This is testing for sig gender differences in controls only for qual traits.
	if (xpar::af_sex_test && !par::af_read){
		if (par::qt) 
			error("Frequency difference test between sexes is for qualitative phenotype only. Quantitative phenotype is detected and command --freqdiff-x is not valid." );
		ofstream FRQ;
		string f = par::output_file_name + ".xtest";
		if (par::summ_nonfounders)
			P.printLOG("Writing X allele frequencies differences (all individuals) to [ " + f + " ] \n");
		else
			P.printLOG("Writing X allele frequencies differences (founders-only) to [ " + f + " ] \n");
		FRQ.open(f.c_str(), ifstream::out);
		FRQ.precision(4);
		FRQ << setw(4) << "CHR" << " "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(4) << "A1" << " "
		<< setw(4) << "A2" << " "
		<< setw(8) << "M1" << " "
		<< setw(8) << "M2" << " "
		<< setw(8) << "F1" << " "
		<< setw(8) << "F2" << " "
		<< setw(12) << "Pvalue"
		<< "\n";

		vector<bool>::iterator d = del.begin();
		vector<Locus*>::iterator loc = P.locus.begin();
		vector< string >::iterator a1_it = allele_1.begin();
		vector< string >::iterator a2_it = allele_2.begin();
		vector< double >::iterator freq_f_it = freq_f.begin(); 
		vector< double >::iterator freq_m_it = freq_m.begin(); 
		vector< int >::iterator nm_f_it = nm_f.begin(); 
		vector< int >::iterator nm_m_it = nm_m.begin(); 
		while (loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){ //only if this locus is on chrX
				string a1 = *a1_it;
				string a2 = *a2_it;
				if (a1=="") a1="0";
				if (a2=="") a2="0";
				FRQ << setw(4)  << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp) << (*loc)->name  << " "
				<< setw(4)  << a1  << " "
				<< setw(4)  << a2  << " "
				<< setw(8)  << *freq_m_it  << " "
				<< setw(8)  << *nm_m_it-*freq_m_it  << " "
				<< setw(8)  << *freq_f_it  << " "
				<< setw(8)  << *nm_f_it-*freq_f_it  << " ";
				table_t t;
				double pvalue;
				sizeTable(t,2,2);
				t[0][0] = *freq_m_it;
				t[0][1] = *nm_m_it - *freq_m_it;
				t[1][0] = *freq_f_it;
				t[1][1] = *nm_f_it - *freq_f_it;
				pvalue = fisher(t);
				if ( pvalue > -1 )	FRQ << setw(12) << pvalue << "\n";
				else	FRQ << setw(12) << "NA" << "\n";
				if (pvalue <= xpar::af_sex_test_limit && pvalue > -1){
					exc_sexdiff++;	
					*d = true;
				}
			}
			loc ++;
			d ++;
			a1_it ++;
			a2_it ++;
			nm_f_it ++;
			nm_m_it ++;
			freq_f_it ++;
			freq_m_it ++;
		}
		FRQ.close();
	}

	///////////////////////////
	// DIFFERENTIAL MISSINGNESS
	// Check for significant differences in genotype missingness rates between 
	// Male and female controls for qualitative traits using Fishers exact test
	if (xpar::am_sex_test){
		// TODO : make this for quantitative traits
		// Should be very similar for each
		if (par::qt) error("Missingness difference test between sexes is for qualitative phenotype only. Quantitative phenotype is detected and command --missdiff-x is not valid." );
		ofstream MIS;
		string m = par::output_file_name + ".xmisstest";
		if (par::summ_nonfounders)
			P.printLOG("Writing X missingness differences (all individuals) to [ " + m + " ] \n");
		else
			P.printLOG("Writing X missingness differences (founders-only) to [ " + m + " ] \n");
		MIS.open(m.c_str(), ifstream::out);
		MIS.precision(4);
		MIS << setw(4) << "CHR" << " "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(4) << "A1" << " "
		<< setw(4) << "A2" << " "
		<< setw(8) << "M1" << " "
		<< setw(8) << "M2" << " "
		<< setw(8) << "F1" << " "
		<< setw(8) << "F2" << " "
		<< setw(12) << "Pvalue"
		<< "\n";

		// TODO : Implement differential missingness filtration
		vector<bool>::iterator d = del.begin();
		vector<Locus*>::iterator loc = P.locus.begin();
		vector< string >::iterator a1_it = allele_1.begin();
		vector< string >::iterator a2_it = allele_2.begin();
		vector< double >::iterator miss_f_it = miss_f.begin(); 
		vector< double >::iterator miss_m_it = miss_m.begin(); 
		vector< int >::iterator nm_f_it = nm_f.begin(); 
		vector< int >::iterator nm_m_it = nm_m.begin(); 
		while (loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){ //only if this locus is on chrX
				string a1 = *a1_it;
				string a2 = *a2_it;
				if (a1=="") a1="0";
				if (a2=="") a2="0";

				MIS << setw(4)             << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp)    << (*loc)->name  << " "
				<< setw(4)  << a1          << " "
				<< setw(4)  << a2          << " "
				<< setw(8)  << *miss_m_it  << " "
				<< setw(8)  << *nm_m_it   << " "
				<< setw(8)  << *miss_f_it  << " "
				<< setw(8)  << *nm_f_it << " ";
				table_t t;
				double pvalue;
				sizeTable(t,2,2);
				t[0][0] = *miss_m_it;
				t[0][1] = *nm_m_it;
				t[1][0] = *miss_f_it;
				t[1][1] = *nm_f_it;
				pvalue = fisher(t);
				if ( pvalue > -1 ) MIS << setw(12) << pvalue << "\n";
				else	MIS << setw(12) << "NA" << "\n";
				if (pvalue <= xpar::am_sex_test_limit && pvalue > -1){
					exc_sexdiff++;	
					*d = true;
				}
			}
			loc ++;
			d ++;
			a1_it ++;
			a2_it ++;
			nm_f_it ++;
			nm_m_it ++;
			miss_f_it ++;
			miss_m_it ++;
		}

		MIS.close();

	}

	/////////////////////////
	// Write HWE statistics

	if (par::HWD_test || par::HWD_report){
		ofstream HWD;
		if (par::HWD_report){
			if (par::summ_nonfounders)
				P.printLOG("Writing Hardy-Weinberg tests (all individuals) to [ " + par::output_file_name + ".hwe ] \n");
			else
				P.printLOG("Writing Hardy-Weinberg tests (founders-only) to [ " + par::output_file_name + ".hwe ] \n");
			string f = par::output_file_name + ".hwe";
			HWD.open(f.c_str(), ifstream::out);
			HWD.precision(4);
			HWD << setw(4) << "CHR" << " "
			<< setw(par::pp_maxsnp) << "SNP" << " "
			<< setw(8) << "TEST" << " "
			<< setw(4) << "A1" << " "
			<< setw(4) << "A2" << " "
			<< setw(20) << "GENO" << " "
			<< setw(8) << "O(HET)" << " "
			<< setw(8) << "E(HET)" << " "
			<< setw(12) << "P" << " "
			<< "\n";
		}

		int cnt=0, cnt_a=0, cnt_u=0;
		////////////////////////
		// Consider each locus

		vector<bool>::iterator d = del.begin();
		vector<Locus*>::iterator loc = P.locus.begin();

		for ( int l = 0 ; l < P.locus.size() ; l++ ){
			// Compute p-values for HWE test in cases, controls & all
			// Only consider founders
			int a11, a12, a22;
			int u11, u12, u22;
			int b11, b12, b22;
			a11=a12=a22=0;
			u11=u12=u22=0;
			b11=b12=b22=0;

			bool X = false, haploid = false;
			if (par::chr_sex[(*loc)->chr])	X=true;
			else if (par::chr_haploid[(*loc)->chr])	haploid=true;
			///////////////////////////////////////////////
			// Iterate over each individual, founders only

			for ( int i = 0 ; i < P.sample.size() ; i++ ){
				Individual * person = P.sample[i];
				///////////////////////////////////////////////
				// Only consider founders, & diploid genotypes
				if ( par::summ_nonfounders || person->founder )
					if ( ! ( haploid || ( X && person->sex ) ) ){
						bool s1 = P.SNP[l]->one[i];
						bool s2 = P.SNP[l]->two[i];
						// Consider everybody, irrespective of phenotype
						// (QT, C/C or missing)
						if (!s1){
							if (!s2) b11++;   //   00 = hom(11)
							else b12++;       //   01 = het(12)
						}
						else if ( s2 )	b22++; // 11 = hom(22)
						if (par::bt){  // for binary trait, separately for cases/controls
							if (person->phenotype == 1){
								if (!s1){
									if (!s2) u11++;   //   00 = hom(11)
									else u12++;       //   01 = het(12)
								}
								else if ( s2 ) u22++; //   11 = hom(22)
							}
							else if (person->phenotype == 2){
								if (!s1){
									if (!s2) a11++;   //   00 = hom(11)
									else a12++;         //   01 = het(12)
								}
								else if ( s2 ) a22++; //   11 = hom(22)
							}
						}
					}
			// Next individual
			}

			// Allele frequencies
			double afreq = 0, ufreq = 0, freq = 0;
			bool include_cases = true;
			bool include_controls = true;
			if (par::qt)	freq = ( b11 + (double)b12/2.0 ) / (double)( b11+b12+b22 );
			else{
				afreq = ( a11 + (double)a12/2.0 ) / (double)( a11+a12+a22 );
				ufreq = ( u11 + (double)u12/2.0 ) / (double)( u11+u12+u22 );
				freq =  ( b11 + (double)b12/2.0 ) / (double)( b11+b12+b22 );
				if ( a11+a12+a22 == 0 ) include_cases = false;
				if ( u11+u12+u22 == 0 ) include_controls = false;
			}
			if (par::qt){
				double p;
				if (par::HWD_standard){
				double tot = b11 + b12 + b22;
				double exp_11 = freq * freq * tot;
				double exp_12 = 2 * freq * (1-freq) * tot;
				double exp_22 = (1-freq) * (1-freq) * tot;
				double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
							+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
							+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;
				p = chiprobP(chisq,1);
				}
				else	p = SNPHWE( b12, b11, b22 );
				if (par::HWD_report){
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name << " "
						<< setw(8) << "ALL(QT)" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20) << (int2str(b11)+"/"+int2str(b12)+"/"+int2str(b22)) << " "
						<< setw(8) << (double)b12/(double)(b11+b12+b22) << " "
						<< setw(8) << 2 * freq * (1-freq)  << " ";
					if ( realnum(p) )	HWD << setw(12) << p << "\n";
					else	HWD << setw(12) << "NA" << "\n";
				}
				if ( p <= par::HWD_limit && p > -1 ){
					cnt++;
					*d = true;
				}
			}
			else{
				// For case/control data
				double p, p_a, p_u;
				if (par::HWD_standard){
					double exp_a11 = afreq * afreq * (a11+a12+a22);
					double exp_a12 = 2 * afreq * (1-afreq) * (a11+a12+a22);
					double exp_a22 = (1-afreq) * (1-afreq) * (a11+a12+a22);
					double exp_u11 = ufreq * ufreq * (u11+u12+u22);
					double exp_u12 = 2 * ufreq * (1-ufreq) * (u11+u12+u22);
					double exp_u22 = (1-ufreq) * (1-ufreq) * (u11+u12+u22);
					double exp_11 = freq * freq * (b11+b12+b22);
					double exp_12 = 2 * freq * (1-freq) * (b11+b12+b22);
					double exp_22 = (1-freq) * (1-freq) * (b11+b12+b22);
					double chisq_a = ( (a11-exp_a11)*(a11-exp_a11) ) / exp_a11
									+ ( (a12-exp_a12)*(a12-exp_a12) ) / exp_a12
									+ ( (a22-exp_a22)*(a22-exp_a22) ) / exp_a22 ;
					double chisq_u = ( (u11-exp_u11)*(u11-exp_u11) ) / exp_u11
									+ ( (u12-exp_u12)*(u12-exp_u12) ) / exp_u12
									+ ( (u22-exp_u22)*(u22-exp_u22) ) / exp_u22 ;
					double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
									+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
									+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;
					p = chiprobP(chisq,1);
					p_a = chiprobP(chisq_a,1);
					p_u = chiprobP(chisq_u,1);
				}
				else{
					p = SNPHWE( b12, b11, b22 );
					p_a = SNPHWE( a12, a11, a22 );
					p_u = SNPHWE( u12, u11, u22 );
				}
				if (par::HWD_report){
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name  << " "
						<< setw(8) << "ALL" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20)
						<< int2str(b11)+"/"+int2str(b12)+"/"+int2str(b22) << " "
						<< setw(8) << (double)b12/(double)(b11+b12+b22) << " "
						<< setw(8) << 2 * freq * (1-freq)  << " ";
					if ( p > -1 )	HWD << setw(12) << p  << "\n";
					else	HWD << setw(12) << "NA"  << "\n";
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name  << " "
						<< setw(8) << "AFF" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20)
						<< int2str(a11)+"/"+int2str(a12)+"/"+int2str(a22) << " "
						<< setw(8) << (double)a12/(double)(a11+a12+a22) << " "
						<< setw(8) << 2 * afreq * (1-afreq)  << " ";
					if (include_cases && p_a > -1 )	HWD << setw(12) << p_a  << "\n";
					else	HWD << setw(12) << "NA" << "\n";
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name  << " "
						<< setw(8) << "UNAFF" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20)
						<< int2str(u11)+"/"+int2str(u12)+"/"+int2str(u22) << " "
						<< setw(8) << (double)u12/(double)(u11+u12+u22) << " "
						<< setw(8) << 2 * ufreq * (1-ufreq)  << " ";
					if (include_controls && p_u > -1 )	HWD << setw(12) << p_u  << "\n";
					else	HWD << setw(12) << "NA" << "\n";
				}
				// Increase counts: in cases
				if ( include_cases && p_a < par::HWD_limit && p_a > -1 ) cnt_a++;
				// Controls (and, if possible, exclude on this value)
				if ( include_controls && p_u < par::HWD_limit && p_u > -1 ){
					cnt_u++;
					if ( ! par::HWD_filter_on_all ){
						*d = true;
						cnt++;
					}
				}
				// In total sample, and if needed, exclude here
				if ( p < par::HWD_limit && p>-1 ){
					if ( par::HWD_filter_on_all || ! include_controls ){
						*d = true;
						cnt++;
					}
				}
	    	}
	    	// next locus
			++loc;
			++d;
		}
		// Finish the report...
		if (par::HWD_report)	HWD.close();
		// ...or finish pruning
		P.printLOG( int2str(cnt) + " markers to be excluded based on HWE test ( p <= " + dbl2str(par::HWD_limit) + " )\n");
		if (par::bt){
			P.printLOG("\t" + int2str(cnt_a) + " markers failed HWE test in cases\n");
			P.printLOG("\t" + int2str(cnt_u) + " markers failed HWE test in controls\n");
		}
	}
	
	///////////////////////////////////////////////////
	// Summary statistics for genotyping/missing rates
	if (par::report_missing){
		///////////////////////////////////////////
		// Report by genotyping rate by individual
		// possibly allowing for obligatory missingness
		P.printLOG( "Writing individual missingness information to [ " + par::output_file_name + ".imiss ] \n");
		ofstream MIS;
		string f = par::output_file_name + ".imiss";
		MIS.open(f.c_str(), ifstream::out);
		MIS.precision(4);
		MIS << setw(par::pp_maxfid) << "FID" << " "
			<< setw(par::pp_maxiid) << "IID" << " "
			<< setw(10) << "MISS_PHENO" << " "
			<< setw(8) << "N_MISS" << " ";
		MIS << setw(8) << "N_GENO" << " ";
		MIS << setw(8) << "F_MISS" << "\n";
		for (int i=0; i<P.n; i++){
			MIS << setw(par::pp_maxfid) << P.sample[i]->fid << " "
				<< setw(par::pp_maxiid) << P.sample[i]->iid << " ";
			if (P.sample[i]->missing) MIS << setw(10) << "Y" << " ";
			else MIS << setw(10) << "N" << " " ;
			// Sum missingness over all SNPs
			int m=0;       // Missing SNPs
			int nsnps=0;   // All non-obligatory missing SNPs
			bool female = ! (P.sample[i]->sex);
			if ( ! par::oblig_missing ){
				for (int l=0; l<P.locus.size();l++){
					// Skip female Y chromosomes
					if ( female && par::chr_Y[P.locus[l]->chr] )	continue;
					if ( MISSING(i,l) ) ++m;
					++nsnps;
				}
			}
			else{ // ... allow oblig missing values
				for (int l=0; l<P.locus.size();l++){
					// Skip female Y chromosomes
					if ( female && par::chr_Y[P.locus[l]->chr] )	continue;
					if ( ! (P.obligMissing(i,l)) ){
						if ( MISSING(i,l) )	++m;
						++nsnps;
					}
				}
			}
			MIS << setw(8) << m << " ";
			MIS << setw(8) << nsnps << " ";
			MIS << setw(8) << (double)m/(double)nsnps << "\n";
		}
		MIS.close();

		///////////////////////////////////////////
		// Report by genotyping rate by locus
		// possibly allowing for sample strata
		// possibly allowing for obligatory missingness
		P.printLOG("Writing locus missingness information to [ " + par::output_file_name +".lmiss ] \n");
		f = par::output_file_name + ".lmiss";
		MIS.open(f.c_str(), ifstream::out);
		MIS.clear();
		MIS.precision(4);

		MIS << setw(4) << "CHR" << " " << setw(par::pp_maxsnp) << "SNP" << " ";
		if (par::include_cluster_from_file)	MIS << setw(10) << "CLST" << " ";
		MIS << setw(8) << "N_MISS" << " ";
		MIS << setw(8) << "N_GENO" << " ";
		if (par::include_cluster_from_file)	MIS << setw(8) << "N_CLST" << " ";
		MIS << setw(8) << "F_MISS" << "\n";
		for (int l=0; l<P.locus.size(); l++){
			Locus * loc = P.locus[l];
			bool chrY = par::chr_Y[P.locus[l]->chr];
			// nk==1 for basic missingness (i.e. not stratified by cluster)
			for (int k=0; k<P.nk; k++){
				MIS << setw(4) << loc->chr << " " << setw(par::pp_maxsnp) << loc->name << " ";
				if (par::include_cluster_from_file)	MIS << setw(10) << P.kname[k] << " ";
				int m=0;     // Number of missing genotypes
				int c=0;     // Number of people in cluster
				int nsnps=0; // Number of actual genotypes in cluster
				for ( int i=0; i<P.sample.size(); i++){
					// Skip female Y chromosome calls
					if ( chrY && ! P.sample[i]->sex )	continue;
					if (par::include_cluster_from_file){
						if ( P.sample[i]->sol == k ){
							if ( ( ! par::oblig_missing ) || ( ! (P.obligMissing(i,l)) ) ){
								if ( MISSING(i,l) ) ++m;
								++nsnps;
							}
							++c;
						}
					}
					else{ // ... ignore cluster strata
						if ( ( ! par::oblig_missing ) || ( ! (P.obligMissing(i,l)) ) ){
							if ( MISSING(i,l) ) ++m;
							++nsnps;
						}
					}
				// Next individual
				}
				MIS << setw(8) << m << " ";
				if (par::include_cluster_from_file)	MIS << setw(8) << c << " ";
				MIS << setw(8) << nsnps << " ";
				MIS << setw(8) << (double)m / (double)nsnps << "\n";
			}
		// Next SNP
		}
		MIS.close();
	}

	/////////////////////////////////
	// Remove rare SNPs

	loc = P.locus.begin();
	d = del.begin();
	while ( loc != P.locus.end() ){
		// Note epsilon correction for MAF, due to floating point
		// issues: only apply to the lower MAF range
		if ( (*loc)->freq < 0 || (*loc)->freq + par::epsilon < par::min_af || (*loc)->freq > par::max_af ){
			*d = true;
			exc_maf++;
		}
		d++;
		loc++;
	}

	/////////////////////////////////////////
	// Remove SNPs based on thresholds

	if ( P.locus.size() > 0 ) P.printLOG("Total genotyping rate in remaining individuals is " + dbl2str(total_genotyping/(double)P.locus.size())+"\n");
	P.printLOG(int2str(exc_miss) + " SNPs failed missingness test ( GENO > " + dbl2str(par::MAX_GENO_MISSING)+" )\n");
	P.printLOG(int2str(exc_maf)+" SNPs failed frequency test ( MAF < "+dbl2str(par::min_af));
	if (par::max_af < 0.5 ) P.printLOG(" or MAF > " + dbl2str(par::max_af));
  	P.printLOG(" )\n");
  	if (xpar::af_sex_test) //added by Feng Gao & Diana Chang. Output this sentence only if the test is specified.
		P.printLOG(int2str(exc_sexdiff)+" SNPs failed sex frequency difference test ( p-value < " + dbl2str(xpar::af_sex_test_limit)+" )\n");
  	
	int tmp = P.deleteSNPs(del);
	//////////////////////////////////////////
	// Need to make back to individual major?

	if ( ! original_SNP_major ) P.SNP2Ind();
	return;
}
