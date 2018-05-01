
///////////////////////////////////////////////////////////////////////////
//       Added by Arjun Biddanda. Modified by Liang Zhang.               //
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

// Note that if you only have one heterozygous female you 
// cannot do this test
double calcVar(vector_t& z, double mu){
	double n = z.size();
	double var = 0;
	for (vector<double>::iterator it = z.begin(); it != z.end(); ++it){
		var += square((*it - mu));
	}
	var /= (n - 1.);
	return var;
}

double calcT(double z1, double z02, double s1, double s02, double n1, double n02){
	double t = (abs(z1 - z02)) / sqrt(s1/n1 + s02/n02);
	return t;
}

double calcDF(double s1, double s02, double n1, double n02){
	double num = (s1 / n1) + (s02 / n02);
	num = square(num);
	double denom = square((s1 / n1)) / (n1 - 1.0) + square((s02 / n02)) / (n02 - 1.0);
	double df = num / denom;
	return df;
}

matrix_t getPhenotypes(Plink& P){ //basically the same as Plink::glmAssoc(), except that this function returns phenotype values
	matrix_t phenotypes;
	vector<Locus *>::iterator loc = P.locus.begin();
	for ( loc ; loc < P.locus.end(); loc++ ){
		vector_t phenotemp(P.n, -1);
		for (int i = 0 ; i < P.n ; i++ ){
			phenotemp[i] = P.sample[i]->pperson->phenotype;
		}
		phenotypes.push_back(phenotemp);
	}
	return phenotypes;
}

double calcWeight(vector_t& y, double mu){
	double n = y.size();
	double weight = 0;
	double var = 0;
	for (vector<double>::iterator iter = y.begin(); iter != y.end(); ++iter) {
		var += square((*iter - mu));
	}
	var /= (n-1); //var = empirical sample variance edited by Liang
	if ( var == 0 ) { 
		weight = DBL_EPSILON;
	}
	else {
		weight = 1/var; 
	}

	return weight;
}

// TODO : aab227 - try to make this more efficient?
vector_t calcWeighting(Plink& P, int SNP, boolmatrix_t ones, boolmatrix_t twos, matrix_t phenotypes) {
    vector_t weights;
	vector_t y00, y01, y11;

    for(int i = 0; i < P.n; i++) {
      bool allele1 = ones[SNP][i];
      bool allele2 = twos[SNP][i];
      Individual* person = P.sample[i];

      if ( !(person->missing) && !(person->sex) ){
        if ( !allele1 ) {
          if ( !allele2 ) {
            y00.push_back(phenotypes[SNP][i]);
          }
          else {
            y01.push_back(phenotypes[SNP][i]);
          }
        }
        else if ( allele2 ) {
          y11.push_back(phenotypes[SNP][i]);
        }
      }
    }

    double y00_mean = calcMean(y00);
    double y01_mean = calcMean(y01);
    double y11_mean = calcMean(y11);
		double w00 = calcWeight(y00, y00_mean);
		double w01 = calcWeight(y01, y01_mean);
		double w11 = calcWeight(y11, y11_mean);

    for(int i = 0; i < P.n; i++) {
    	bool allele1 = ones[SNP][i];
    	bool allele2 = twos[SNP][i];
    	Individual* person = P.sample[i];

        if ( !(person->missing) && !(person->sex) ){
            if ( !allele1 ) {
                if ( !allele2 ) {
                    weights.push_back(w00);
                }
                else {
                    weights.push_back(w01);
                }
            }
            else if ( allele2 ) {
                weights.push_back(w11);
            }
        }
    }

    return weights;
}

// basically the same as Plink::glmAssoc(), except that this function returns residuals
matrix_t xglmAssocResidualsVarHet(Plink& P, boolmatrix_t ones, boolmatrix_t twos, matrix_t phenotypes, int flag){ 
	matrix_t residuals; //stores the residuals 
	Plink* Q = &P;
	P.SNP2Ind();
	vector<Locus*>::iterator loc = P.locus.begin();
	
	// Make sure to reset missing status for real afterwards..
	vector< bool > missingStatus(P.sample.size());
	vector<Individual*>::iterator person = P.sample.begin();
	vector< bool >::iterator missing_it = missingStatus.begin();
	while(person != P.sample.end()){ //retain a copy for the real missing status
		*missing_it = (*person)->missing;
		person ++;
		missing_it ++;
	}

	// Set all males as missing
	for (int i=0; i<P.n; i++){
		if ( ! P.sample[i]->missing ){
			if ( P.sample[i]->sex )	P.sample[i]->missing = true;
		}
	}

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
				// Build design matrix
				lm->buildDesignMatrix();
				//////////////////////////////
				// Clusters specified?
				if ( par::include_cluster ){
					lm->setCluster();
				}

				//////////////////////////////////////////////////
				// Fit linear or logistic model (Newton-Raphson)
				vector< vector<double> > Xdes;
				int nind = lm->Ysize();
				int np = lm->getNP();

				if ( flag == 0 ){
					// cout << "Running for SNP : " << (*loc)->name << "\n";

					lm->fitLM();
					lm->validParameters();

					Xdes = lm->X;

					vector_t resid_temp(nind, -1);
					if (lm->isValid() && !(lm->getVar()[1] < 1e-20) && realnum(lm->getVar()[1])){ //the model or result is OK
						// i = index in nind , j = index in P.n
						int i = 0;
						int j = 0;
						while (i < nind && j < P.n){	
							if (!(P.sample[j]->missing)){
								double e_i = P.sample[j]->pperson->phenotype;
								for (int p=0; p<np; p++){
									e_i -= (lm->getCoefs())[p] * Xdes[i][p];
								}
								e_i = abs(e_i);
								resid_temp[i] = e_i;
								i++;
							}
							j++;
						}
					}
					residuals.push_back(resid_temp);

					delete lm;
				}

				else if(flag == 1){
					//save original phenotype values to restore later
					vector_t orig_phenotypes(P.n, -1);
					for(int i = 0; i < P.n; i++) {
						orig_phenotypes[i] = P.sample[i]->pperson->phenotype;
					}
					vector_t weights = calcWeighting(P, l, ones, twos, phenotypes);
					for(int a = 0; a < lm->X.size(); a++){
						for(int b = 0; b < lm->X[a].size(); b++){
							if(weights[a] > 0) {
								lm->X[a][b] = lm->X[a][b] * sqrt(weights[a]);
							}
						}
					}

					// TODO : lb598 - lm is not fit here with new design parameters... why?
					// TODO : lb598 - should weight phenotypes before lm->fitLM();

					int i = 0;
					int j = 0;
					while (i < nind && j < P.n){
						if (!(P.sample[j]->missing)){
                                                       if (weights[i] > 0){
                                                             P.sample[j]->pperson->phenotype = orig_phenotypes[j] * sqrt(weights[i]);
							}
							i++;
						}
						j++;
					}
					lm->setDependent();
					lm->fitLM();
					lm->validParameters();
					lm->isValid();
					vector_t pvals = lm->getPVals();
					//note here residuals actually are the p values for the weighted linear regressions
					residuals.push_back(pvals); 
					//restore phenotype values
					for(int i = 0; i < P.n; i++) {
						P.sample[i]->pperson->phenotype = orig_phenotypes[i];
					}

					delete lm;
				}

			}
			loc ++;
			l ++;
		}
	}

	person = P.sample.begin();
	missing_it = missingStatus.begin();
	while(person != P.sample.end()){
		(*person)->missing = *missing_it;
		person ++;
		missing_it ++;
	}//recover P missing statuses

	par::assoc_glm_without_main_snp = orig_without_main_snp;
	par::standard_beta = orig_standard_beta;

	P.Ind2SNP(); // Lauren 3/4/18: restore integrity of P

	return residuals;
}

void xVarianceHeterogeneity(Plink& P, vector<double>& pval_results) {
	if (!par::qt) {
		if (xpar::var_het == true) error("--var-het for quantitative traits only");
		if (xpar::var_het_weight == true) error("--var-het-weight for quantitative traits only");
		if (xpar::var_het_comb == true) error("--var-het-comb for quantitative traits only");
	}
	string f = par::output_file_name;
	//temp file for gene based testing
	string f1 = f + "_varhet_temp.txt";
	string f2 = f + "_varhetweight_temp.txt";
	string f3 = f + "_varhetcomb_temp.txt";
	f += ".xvarhet";
	P.printLOG("Writing variance heterogeneity test results to [ " + f + " ] \n");
	ofstream ASC;
	ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " "
	<< setw(par::pp_maxsnp) << "SNP" << " "
	<< setw(10) << "BP" << " "
	<< setw(4) << "A1" << " ";
	ofstream ASC1;
	ofstream ASC2;
	ofstream ASC3;

	if(xpar::var_het) {
		ASC << setw(10) << "T-STAT" << " "
		<< setw(10) << "DF" << " "
		<< setw(10) << "PVAL" << " ";
		ASC1.open(f1.c_str(), ios::out); //file with snp names [column1] and p values [column2]
	}
	if(xpar::var_het_weight){
		ASC << setw(16) << "T-STAT_WEIGHT" << " "
		<< setw(10) << "DF_WEIGHT" << " "
		<< setw(14) << "PVAL_WEIGHT" << " ";
		ASC2.open(f2.c_str(), ios::out); //file for var_het_weight_gene_based_test

	}
	if(xpar::var_het_comb){
		ASC << setw(10) << "Z-STAT" << " "
		<< setw(10) << "PVAL_COMB" << " ";
		ASC3.open(f3.c_str(), ios::out); //file for var_het_comb_gene_based_test
	}
	ASC << "\n";
	vector<Locus *>::iterator loc1 = P.locus.begin();
	boolmatrix_t ones;
	boolmatrix_t twos;
	for (int l=0; l < P.locus.size(); l++){
		if (par::chr_sex[(*loc1)->chr]){
			ones.push_back((P.SNP[l]->one));
			twos.push_back((P.SNP[l]->two));
		}
		else
			P.printLOG("At least one of your loci is not on X chromosome!");
		++loc1;
	}

	// Get all of the regression residuals for this PLINK fileset
	matrix_t phenotypes = getPhenotypes(P);
	matrix_t residuals, weightedResiduals; 
	//note the weightedResiduals actually means the p values the t statistic for each coef. 
	if(xpar::var_het || xpar::var_het_comb){
		residuals = xglmAssocResidualsVarHet(P, ones, twos, phenotypes, 0);
	}

	if(xpar::var_het_weight || xpar::var_het_comb){
		weightedResiduals = xglmAssocResidualsVarHet(P, ones, twos, phenotypes, 1);
	}

	int xsnpindex = 0;
	vector<Locus *>::iterator loc = P.locus.begin();
	for (int l = 0; l < P.locus.size(); l++){
		vector_t z1, z02, y00, y01, y11, z1_weight, z02_weight;
		double z1_mean, z02_mean, s1, s02, n1, n02, t, df, pval;
		double pval_weight;
		double zscore1, zscore2, zstat, pval_comb;
		if (par::chr_sex[(*loc)->chr]){
			// Need to write output here
		        // CHR, SNP, BP, A1
			ASC << setw(4)  << (*loc)->chr << " "
			<< setw(par::pp_maxsnp) << (*loc)->name << " "
			<< setw(10) << (*loc)->bp << " "
			<< setw(4)  << (*loc)->allele1 << " ";
		}

		if(xpar::var_het || xpar::var_het_comb){

			int j = 0;
			int i = 0;
			//for individual, construct vectors of its phenotype values across all SNPs
			if (residuals[xsnpindex].size() > 0){
				while (j < P.n) {
					Individual * person = P.sample[j];
					if (!(person->missing) && !(person->sex)) { // not male && not missing
						bool s1 = ones[l][j];
						bool s2 = twos[l][j];

						if (residuals[xsnpindex][i] > 0){
							if (!s1){ 
								if (!s2){ //00 == hom(11)
									z02.push_back(residuals[xsnpindex][i]);
								}
								else{ //10 || 01 == het
									z1.push_back(residuals[xsnpindex][i]);
								}
							}	
							else if (s2) {
								z02.push_back(residuals[xsnpindex][i]);
							}
						}
						i++;
					}
					j++;
				}
			}
			n1 = z1.size();
			n02 = z02.size();
			// Set guard so we do not accidentally compute statistic...
			if (n1 > 1 && n02 > 2){
				//Calculate the actual statistic
				z1_mean = calcMean(z1);
				z02_mean = calcMean(z02);
				s1 = calcVar(z1, z1_mean);
				s02 = calcVar(z02, z02_mean);
				t = calcT(z1_mean, z02_mean, s1, s02, n1, n02);
				df = calcDF(s1, s02, n1, n02);
				pval = calc_tprob(t, df);
			}
			else{
				P.printLOG("At locus "+(*loc)->name+" The number of heterozygous female is less than 1 or the number of homozogous female or male is less than 2.\n");
			}
		}

		if(xpar::var_het_weight || xpar::var_het_comb){
			pval_weight = weightedResiduals[l][0]; //it is the p value for beta_1
		}
		// Amended 
		if(xpar::var_het_comb && n1 > 1 && n02 > 2){
			zscore1 = fabs(ltqnorm(pval * 0.5));
			zscore2 = fabs(ltqnorm(pval_weight * 0.5));
			zstat = (zscore1 + zscore2) / sqrt(2);
			pval_comb = 2 * (1 - normdist(fabs(zstat)));
		}
		if(xpar::var_het) {
			if (n1 > 1 && n02 > 2){
				ASC << setw(10) << t << " "
				<< setw(10) << df << " " 
				<< setw(10) << pval << " ";
				ASC1 << setw(par::pp_maxsnp) << (*loc)->name << " "<< setw(10) << pval<<"\n";
				pval_results.push_back(pval);
			}
			else{
				ASC << setw(10) << "N/A" << " "
				<< setw(10) << "N/A" << " " 
				<< setw(10) << "N/A" << " ";
			}
		}
		// Also need to check if the pvalue is from a valid linear model...
		if(xpar::var_het_weight) {
				// TODO : possibly have something other than N/A there...
				ASC << setw(16) << "N/A" << " "
				<< setw(10) << "N/A" << " "
				<< setw(14) << pval_weight << " ";
				ASC2 << setw(par::pp_maxsnp) << (*loc)->name << " "<< setw(14) << pval_weight<<"\n";
		}
		if(xpar::var_het_comb) {
			if(n1 > 1 && n02 > 2) {
				ASC << setw(10) << zstat << " "
				<< setw(10) << pval_comb << " ";
				ASC3 << setw(par::pp_maxsnp) << (*loc)->name << " "<< setw(10) << pval_comb<<"\n";
			}
			else { 
				ASC << setw(10) << "N/A" << " "
				<< setw(10) << "N/A" << " ";
			}
		}
		ASC << "\n";
		xsnpindex++;
		++loc;
	}

}

