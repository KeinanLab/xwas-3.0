
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

void xAssocFunctions(Plink& P){ //only one function now. May add others in the future.
	vector<double> pval;
	if (xpar::var_het || xpar::var_het_weight || xpar::var_het_comb) xVarianceHeterogeneity(P,pval);
	if (xpar::var_het && par::use_GC) {
	  string f = par::output_file_name + ".xvarhet";
	  xMultComp(P,pval,f);
	}
	if (xpar::strat_sex || xpar::sex_diff)	xSexRelatedTests(P);
}

// NOTE: this function needs P to be in individual major format
// Before using, call P.SNP2Ind()
// basically the same as Plink::glmAssoc(), except that this function returns p-values
matrix_t xglmAssoc(Plink& P) { 

	matrix_t pval; // stores both pval and beta
	Plink* Q = &P;

	vector<Locus*>::iterator loc = P.locus.begin();
	/////////////////////////////
	// Determine sex distribution
	int nmales = 0, nfemales = 0;
	for (int i=0; i<P.n; i++)
		if ( ! P.sample[i]->missing ){
			if ( P.sample[i]->sex )	nmales++;
			else	nfemales++;
		}
	bool variationInSex = nmales > 0 && nfemales > 0;
	bool orig_without_main_snp = par::assoc_glm_without_main_snp;
	bool orig_standard_beta = par::standard_beta;
	par::assoc_glm_without_main_snp = false;
	par::standard_beta = false;
	if (par::xchr_model != 0){ //only if chrX is considered
		int l = 0;
		while ( loc != P.locus.end() ){
     		if (par::chr_sex[(*loc)->chr] || par::chr_haploid[(*loc)->chr]){ //only if chrX SNP or haplotype to save running time
				bool X=true;
				bool automaticSex=false;
				Model * lm;
				if (par::bt){
					LogisticModel * m = new LogisticModel(Q);
					lm = m;
				}
				else{
					LinearModel * m = new LinearModel(Q);
					lm = m;
				}
				lm->setMissing();
				if ( par::glm_dominant )	lm->setDominant();
				else if ( par::glm_recessive || par::twoDFmodel_hethom )
					lm->setRecessive();
				string mainEffect = "";
				bool genotypic = false;
				if ( ! par::assoc_glm_without_main_snp ) {//need to assume that this test is for all loci?
					genotypic = par::chr_haploid[(*loc)->chr] ? false : par::twoDFmodel ;
					if ( par::glm_recessive )	mainEffect = "REC";
					else if ( par::glm_dominant )	mainEffect = "DOM";
					else if ( par::twoDFmodel_hethom )	mainEffect = "HOM";
					else	mainEffect = "ADD";	  
					lm->addAdditiveSNP(l); 
					lm->label.push_back(mainEffect);
					if ( genotypic ){
						lm->addDominanceSNP(l);	      
						if ( par::twoDFmodel_hethom )	lm->label.push_back("HET");
						else	lm->label.push_back("DOMDEV");
					}
				}
				//////////////////////////////////////////////////////////
				// Haplotype test
				if ( par::chap_test ){
					for (int h=1; h < P.whap->current->group.size(); h++){
						lm->addHaplotypeDosage( P.whap->current->group[h] );	    
						lm->label.push_back( "WHAP"+int2str(h+1) );
					}
				}     
				if ( par::proxy_glm ){
					set<int> t1 = P.haplo->makeSetFromMap(P.haplo->testSet);
					lm->addHaplotypeDosage( t1 );
					lm->label.push_back( "PROXY" );	    
				}
				if ( par::test_hap_GLM ){
					set<int>::iterator i = P.haplo->sets.begin();
					while ( i != P.haplo->sets.end() ){
						set<int> t;
						t.insert(*i);
						lm->addHaplotypeDosage( t );
						lm->label.push_back( P.haplo->haplotypeName( *i ) );
						++i;
					}
				}
				//////////////////////////////////////////////////////////
				// Conditioning SNPs?      
				if (par::conditioning_snps){
					if ( par::chap_test ){
						for (int c=0; c<P.conditioner.size(); c++){
							if ( P.whap->current->masked_conditioning_snps[c] ){
								lm->addAdditiveSNP(P.conditioner[c]); 
								lm->label.push_back(P.locus[P.conditioner[c]]->name);
							}
						}
					}
					else{
						for (int c=0; c<P.conditioner.size(); c++){
							lm->addAdditiveSNP(P.conditioner[c]); 
							lm->label.push_back(P.locus[P.conditioner[c]]->name);
						}
					}
				}
				//////////////////////////////////////////////////////////      
				// Sex-covariate (necessary for X chromosome models, unless
				// explicitly told otherwise). For male/female-only tests,
				// clearly there is no variation in sex
				if ( ( par::glm_sex_effect || ( X && !par::glm_no_auto_sex_effect ) ) && variationInSex ){
					automaticSex = true;
					lm->addSexEffect();
					lm->label.push_back("SEX");	  
				}
				//////////////////////////////////////////////////////////
				// Covariates?
				if (par::clist){
					for (int c=0; c<par::clist_number; c++){
						lm->addCovariate(c);
						lm->label.push_back(P.clistname[c]);
					}
				}
				//////////////////////////////////////////////////////////
				// Interactions
				// addInteraction() takes parameter numbers
				// i.e. not covariate codes      
				// 0 intercept
				// 1 {A}
				//   {D}
				//   {conditioning SNPs}
				//   {sex efffect}
				//   {covariates}
				// Allow for interactions between conditioning SNPs, sex, covariates, etc
      	
				////////////////////////////////////////
				// Basic SNP x covariate interaction?       
				// Currently -- do not allow interactions if no main effect 
				// SNP -- i.e. we need a recoding of things here.
				if ( par::simple_interaction && ! par::assoc_glm_without_main_snp ){
	  				// A, D and haplotypes by conditioning SNPs, sex, covariates	  
	  				int cindex = 2;
	  				if ( genotypic )	cindex = 3;
	  				for (int c=0; c<P.conditioner.size(); c++){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"xCSNP"+int2str(c+1));
						if ( genotypic ){
							lm->addInteraction(2,cindex);
							if ( par::twoDFmodel_hethom )	lm->label.push_back("HETxCSNP"+int2str(c+1));	  
							else	lm->label.push_back("DOMDEVxCSNP"+int2str(c+1));	  
						}	      
						cindex++;
					}
					if ( automaticSex ){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"xSEX");	  
						if ( genotypic ){
							lm->addInteraction(2,cindex);
							if ( par::twoDFmodel_hethom )	lm->label.push_back("HETxSEX");
							else	lm->label.push_back("DOMDEVxSEX");
						}	      
						cindex++;
					}
					for (int c=0; c<par::clist_number; c++){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"x"+P.clistname[c]);	  	      
						if ( genotypic ){
							lm->addInteraction(2,cindex);		  
							if ( par::twoDFmodel_hethom )		  
								lm->label.push_back("HETx"+P.clistname[c]); 
							else
								lm->label.push_back("DOMDEVx"+P.clistname[c]); 
						}	      
						cindex++;	      
					}
				}

				//////////////////////////////
				// Fancy X chromosome models      
				if ( X && automaticSex && par::xchr_model > 2 ){ //obviously these fancy models are not used with only females or males
					// Interaction between allelic term and sex (i.e. allow scale of male effect to vary)
					int sindex = 2;
					if ( genotypic )	sindex++;
					sindex += P.conditioner.size(); 
					lm->addInteraction(2,sindex);
					lm->label.push_back("XxSEX");	  	  
					// xchr model 3 : test ADD + XxSEX
					// xchr model 4 : test ADD + DOM + XxSEX
				}
				//////////////////////////////
				// Build design matrix
				lm->buildDesignMatrix();
				//////////////////////////////
				// Clusters specified?
				if ( par::include_cluster ){
					lm->setCluster();
				}
				//////////////////////////////////////////////////
				// Fit linear or logistic model (Newton-Raphson)      
				lm->fitLM();
				lm->validParameters();
				vector_t pval_temp(2, -1.);
				if (lm->isValid() && !(lm->getVar()[1] < 1e-20) && realnum(lm->getVar()[1])){ //the model or result is OK
					pval_temp[0] = lm->getCoefs()[1];
					pval_temp[1] = lm->getPVals()[0];
				}
				pval.push_back(pval_temp);
				delete lm; //clear linear model
			}
			loc ++;
			l ++; 
		}
	}
	par::assoc_glm_without_main_snp = orig_without_main_snp;
	par::standard_beta = orig_standard_beta;

	return pval;
}

void xglmAssocSexSep(Plink& P, matrix_t& resMale, matrix_t& resFemale) {

	P.SNP2Ind();

	if (par::bt)
	  P.printLOG("Fitting logistic model to males only \n");
	else
	  P.printLOG("Fitting linear model to males only \n");

	vector<bool> missingStatus(P.sample.size());
	vector<bool>::iterator missing_it = missingStatus.begin();
	vector<Individual*>::iterator person = P.sample.begin();

	while (person != P.sample.end()) { // retain a copy for the real missing status
		*missing_it = (*person)->missing;
		person ++;
		missing_it ++;
	}

	person = P.sample.begin();
	while (person != P.sample.end()) { // Set all males to missing
		if ((*person)->sexcode != "1")
		  (*person)->missing = true;
		person ++;
	}

	resMale = xglmAssoc(P);

	if (par::bt)
	  P.printLOG("Fitting logistic model to females only \n");
	else
	  P.printLOG("Fitting linear model to females only \n");

	person = P.sample.begin();
	missing_it = missingStatus.begin();
	while (person != P.sample.end()) { // Set all females to missing
		if ((*person)->sexcode != "2")
		  (*person)->missing = true;
		else
		  (*person)->missing = *missing_it;
		person ++;
		missing_it ++;
	}

	resFemale = xglmAssoc(P);

	// Reset missing statuses to original
	person = P.sample.begin();
	missing_it = missingStatus.begin();
	while (person != P.sample.end()) { 
		(*person)->missing = *missing_it;
		person ++;
		missing_it ++;
	}

	P.Ind2SNP(); // Added Lauren 2/23/18, this restores integrity of P structure
}

void xSexRelatedTests(Plink& P) {
        int xchr_model_orig = par::xchr_model;

	if (xpar::multi_xchr_model) {
	  P.printLOG("Analyzing data with males coded 0/1\n");
	  par::xchr_model = 1;
	}
	matrix_t resultpv_male, resultpv_female;
	xglmAssocSexSep(P, resultpv_male, resultpv_female);

	string f = par::output_file_name;

	// Added by Lauren, 3/3/18: concurrent testing of males 0/1 and 0/2 coded
	if (xpar::multi_xchr_model) {
	  vector<double> pval1, pval2;
	  matrix_t resultpv_male_02, resultpv_female_02;

	  par::output_file_name = f + ".01";
	  if (xpar::strat_sex) xStratSexAssocTest(P, resultpv_male, resultpv_female, pval1);
	  if (xpar::sex_diff)  xSexDiffTest(P, resultpv_male, resultpv_female);
	  if (xpar::strat_sex && par::use_GC) {
	    if (par::bt) xMultComp(P,pval1,par::output_file_name + ".xstrat.logistic");
	    else xMultComp(P,pval1,par::output_file_name + ".xstrat.linear");
	  }

	  P.printLOG("Analyzing data with males coded 0/2\n");
	  par::xchr_model = 2;
	  xglmAssocSexSep(P, resultpv_male_02, resultpv_female_02);
	  par::output_file_name = f + ".02";

	  if (xpar::strat_sex) xStratSexAssocTest(P, resultpv_male_02, resultpv_female_02, pval2);
	  if (xpar::sex_diff)  xSexDiffTest(P, resultpv_male_02, resultpv_female_02);
	  if (xpar::strat_sex && par::use_GC) {
	    if (par::bt) xMultComp(P,pval2,par::output_file_name + ".xstrat.logistic");
	    else xMultComp(P,pval2,par::output_file_name + ".xstrat.linear");
	  }
	} else {
	  vector<double> pval;
	  if (xpar::strat_sex) xStratSexAssocTest(P, resultpv_male, resultpv_female, pval);
	  if (xpar::sex_diff)  xSexDiffTest(P, resultpv_male, resultpv_female);
	  if (xpar::strat_sex && par::use_GC) {
	    if (par::bt) xMultComp(P,pval,par::output_file_name + ".xstrat.logistic");
	    else xMultComp(P,pval,par::output_file_name + ".xstrat.linear");
	  }
	}

	par::output_file_name = f;
        par::xchr_model = xchr_model_orig;
}
