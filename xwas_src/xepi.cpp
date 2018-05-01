
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

void calcXepistasis(Plink& P){ 

	Plink* Q = &P;
	if( par::SNP_major) P.SNP2Ind();

	////////////////////////////////////////
	// Set up results files
	ofstream XEPI;
	string f = par::output_file_name;
	if(par::qt){
		f += ".xepi.qt"; // qt: use t test
		XEPI.open (f.c_str(), ios::out);
		P.printLOG("Writing epistasis pairwise results to [ " + f + " ]\n");

		XEPI.precision(11);

		XEPI << setw(4) << "CHR1" << " "
		<< setw(par::pp_maxsnp) << "SNP1" << " "
		<< setw(4) << "CHR2" << " "
		<< setw(par::pp_maxsnp) << "SNP2" << " "
		<< setw(12) << "BETA_INT" << " "
		<< setw(12) << "T_statistic" << " "
		<< setw(12) << "P_value" << " "
		<< "\n";
	} else {
		f += ".xepi.cc"; // bt: use Wald test
		XEPI.open (f.c_str(), ios::out);
		P.printLOG("Writing epistasis pairwise results to [" + f + "]\n");

		XEPI.precision(4);

		XEPI << setw(4) << "CHR1" << " "
		<< setw(par::pp_maxsnp) << "SNP1" << " "
		<< setw(4) << "CHR2" << " "
		<< setw(par::pp_maxsnp) << "SNP2" << " "
		<< setw(12) << "OR_INT" << " "
		<< setw(12) << "Chiq_statistic" << " "
		<< setw(18) << "P_value" << " "
		<< "\n";
	}


	/////////////////////////////////////////////////////////////////////
	// xepi1 and xepi2 threshold were given in terms of 0.0001(two-sides)
	
	P.printLOG("Threshold for displaying xepistatis results (--xepi1): p <=" + dbl2str(xpar::xepi_alpha1) + "\n");
	P.printLOG("Threshold for counting xepistasis results (--xepi2): p <=" + dbl2str(xpar::xepi_alpha2) + "\n");
	
	// Regression based test: case/control or quantitative trait
	
	// Take a list of SNPs, or all SNPs (vector<bool> sA) 
	// Test these against either themselves, or all SNPs (vector<bool> sB)

	// A    x    B
	// ALL  x    ALL   skip e1 > e2
	// SET1 x    ALL  
	// SET1 x    SET1  skip e1 > e2
	// SET1 x    SET2	

	bool skip_symm = false; // To avoid redundant calculation

	// Only output epistatic tests that have p < xpar::xepi_alpha1;
	// Do not even attempt to save any epistatic results

	// Also present summary results for each SNP
	// (i.e. average / proportion of significant epistatic tests)
	// at a certain alpha level, xpar::xepi_alpha2

	vector<bool> sA(P.nl_all, false);
	vector<bool> sB(P.nl_all, false);

	// Are we using a set test? If so, construct now
	if(par::set_test){
		if (P.snpset.size() > 2)
			error ("Can only specify one or two SETs when testing for xepistaisis");
		if (P.snpset.size() == 0)
			error ("There are no valid sets specified");
		for (int e=0;e<P.snpset[0].size();e++){
			sA[P.snpset[0][e]] = true;
		}

		if(P.snpset.size()==2) {
			P.printLOG ("SET1 X SET2 xepistasis mode\n");
			for (int e=0; e<P.snpset[1].size();e++){
				
				sB[P.snpset[1][e]] = true;
			}
		}else if (xpar::Set_by_Set){
			P.printLOG ("SET1 X SET1 xepistasis mode\n");
			skip_symm = true;
			for( int e=0;e<P.snpset[0].size();e++){
				sB[P.snpset[0][e]] = true;
			}
		}else {
			P.printLOG ("SET1 X ALL xepistasis mode\n");
			for (int e=0;e<P.nl_all;e++){
				
				sB[e] = true;
			}
		}

	}else{
		P.printLOG ("ALL X ALL xepistasis mode\n");
		skip_symm = true;
		for (int e=0;e<P.nl_all; e++){
			sA[e] = true;
			sB[e] = true;
		}
	}

	// Count how many items in the SET1

	int epc = 0; // for sA
	int epcb = 0; // for sB
	for (vector<bool>::iterator e1 = sA.begin(); e1 != sA.end();e1++)
		if (*e1) epc++;
	
	for (vector<bool>::iterator e2 = sB.begin(); e2 != sB.end();e2++)
		if (*e2) epcb++;

	int epcc = 0;

	// Keep track of how many xepistasis tests actually performed
	long int nepi = 0;

	vector<int> summary_sig(P.nl_all,0);
	vector<int> summary_good(P.nl_all,0);
	vector<double> best_score(P.nl_all,1);
	vector<int> best_partner(P.nl_all);

	//////////////////////////////////////////
	// Begin interating over pairs: SET X SET
	vector<double> H44_value; // a row vector for recoding H44 value 
	matrix_t corr_sparse; 
	sizeMatrix(corr_sparse, epc*epcb,0); // max possible size
	int corr_index = 0; // The actually rows of corr_sparse, equal to the valid pvalue

        vector<bool> miss_ind(P.n,false);
        int n_miss_ind = 0;
	
	for(int e1=0;e1<P.nl_all;e1++){
		if(sA[e1]){
			for(int e2=0;e2<P.nl_all;e2++){
				//////////////////////////////////////////
				// Skip this test under certain conditions

				if (!sB[e2]) continue; // The SNP is not in the test
                                if (e1>=e2 && skip_symm) continue; // We've already performed this test
                                if (e1==e2) continue; // Same SNP 

                                //////////////////////////////////
                                // Perform test of epistasis here
                                Model * lm;

                                if(par::bt) { // bt: logistic regression with Wald test
                                        LogisticModel * m = new LogisticModel(Q);
					lm = m;

					//Set missing data
					lm -> setMissing();

					//Main effect of SNP1
					lm -> addAdditiveSNP(e1);
					lm -> label.push_back("ADD1");

					//Main effect of SNP2
					lm -> addAdditiveSNP(e2);
					lm -> label.push_back("ADD2");

					//Epistasis
					lm -> addInteraction(1,2);
					lm -> label.push_back("XEPI");

					//Add covariates?
					if (par::clist){
						for (int c=0;c<par::clist_number;c++){
							lm -> addCovariate(c);
							lm -> label.push_back(P.clistname[c]);
						}
					}

					//Build design matrix
					lm -> buildDesignMatrix();

					//Fit regression model
					lm -> fitLM();

					// Did model fit OK?
					lm -> validParameters();

					// Obtain estimates and use Wald test

		   			lm->testParameter = 3; // interaction
		  			vector_t b = lm->getCoefs();
		   			double chisq = lm->getStatistic();
		  			double pvalue = chiprobP(chisq,1);
		   			double z = sqrt(chisq);

		   			if(lm -> isValid()){
						// One more valid test performing
						nepi ++;

						// Count as a good result
						summary_good[e1]++;
						if(sA[e2]) summary_good[e2] ++;

						//Do we want to record this as part of the summary for the first set?
						if (pvalue <= xpar::xepi_alpha2) {
							// first variable will always be in A set
							summary_sig[e1]++;

							// but the sechond may also be in A set
							if(sA[e2]) summary_sig[e2]++;
						}
					}

					// Is this result the best score yet for marker in set A
					if (pvalue <= best_score[e1]){
						best_score[e1] = pvalue;
						best_partner[e1] = e2;
					}

					// The second marker might also be in set A
					if(sA[e2]){
						if(pvalue <= best_score[e2]){
							best_score[e2] = pvalue;
							best_partner[e2] = e1;
						}
					}

					// Is this worth displaying?
					if(pvalue <= xpar::xepi_alpha1){
						XEPI << setw(4) << P.locus[e1]->chr << " "
						<< setw(par::pp_maxsnp) << P.locus[e1]->name << " "
						<< setw(4) << P.locus[e2]->chr << " "
						<< setw(par::pp_maxsnp) << P.locus[e2]->name << " ";

						if (lm -> isValid()){
							XEPI << setw(12) << exp(b[3]) << " "
							<< setw(12) << chisq << " "
							<< setw(12) << pvalue << " "
							<< "\n";
						}else {
							XEPI << setw(12) << "NA" << " "
							<< setw(12) << "NA" << " "
							<< setw(12) << "NA" << " "
							<< "\n";
						}
						XEPI.flush();
					}


	      		} else { // qt: linear regression with t test
                                        LinearModel * m = new LinearModel(Q);
					lm = m;

					//Set missing data
					lm -> setMissing();

					//Main effect of SNP1
					lm -> addAdditiveSNP(e1);
					lm -> label.push_back("ADD1");

					//Main effect of SNP2
					lm -> addAdditiveSNP(e2);
					lm -> label.push_back("ADD2");

					//Epistasis
					lm -> addInteraction(1,2);
					lm -> label.push_back("XEPI");

					//Add covariates?
					if (par::clist){
						for (int c=0;c<par::clist_number;c++){
							lm -> addCovariate(c);
							lm -> label.push_back(P.clistname[c]);
						}
					}

					//Build design matrix
					lm -> buildDesignMatrix();

					if(lm -> isValid()) {

						// Fit regression model
						lm -> fitLM();

						// Did model fit OK?
						lm -> validParameters();

						// Obtain estimate and use t-test
						lm -> testParameter = 3; // parameter for interaction term
						vector_t b = lm -> getCoefs();

						// Since lm excluded indivs w missing phenotype, need to find true design matrix that lm used
						matrix_t Xdes = lm -> X;
						int nind = lm -> Ysize(); // # of non-missing indivs in lm
						int np = lm -> getNP(); // # of columes in Xdes
                                                vector<bool> miss = lm -> getMissing(); // which individuals did lm exclude?

						vector_t resid_squar_temp (nind, -1);
						double mean_b_inter = 0.0; // calculate the mean of interaction column
						int i=0;
						int j=0;
						double resid_sum = 0.0;
						while (j < P.n) {
                                                        if (!miss[j]) {
								double e_i = P.sample[j] -> pperson -> phenotype;
								for (int p=0;p<np;p++){
									e_i -= b[p] * Xdes[i][p];
								}
								mean_b_inter += Xdes[i][3];
								e_i = abs(e_i);

								resid_squar_temp[i] = pow(e_i,2);
								resid_sum += resid_squar_temp[i];
								i++;
							}
							j++;
						}

						int df = 0; // df = #sample - #of parameters except bata0 -1
						if(par::clist) {
							df = nind - 3 - par::clist_number -1;
						} else {
							df = nind -3 -1;
						}
						matrix_t Xdes_t;
                                                transposeMatrix(Xdes,Xdes_t);

						matrix_t XX; // dimension is 4+clist_number
						multMatrix(Xdes_t,Xdes,XX);
						bool flag = true;
						bool all_valid = lm ->checkVIF();
						matrix_t H;
						H = svd_inverse(XX, flag); // Return the inverse of a matrix using singular value decomposition

						double delta = resid_sum/df;
						double t_statistic_inter = b[3]/sqrt(delta*H[3][3]);
						double pvalue = calc_tprob(t_statistic_inter,df);

						if(pvalue <= xpar::xepi_alpha1){
							matrix_t temp_a; // 1*np matrix
							matrix_t H4; // 1*4
							sizeMatrix(H4,1,np);
                                                        H4[0]= H[3];
							multMatrix(H4, Xdes_t,temp_a);

                                                        // create correlation matrix, push NAN for individuals lm excluded
                                                        // keep track of which individuals were excluded from at least one SNP
                                                        int i = 0;
                                                        int j = 0;
                                                        while (j < P.n) {
                                                                if (!miss[j]) {
                                                                        corr_sparse[corr_index].push_back(temp_a[0][i]/sqrt(H[3][3]));
                                                                        i++;
                                                                } else {
                                                                        corr_sparse[corr_index].push_back(NAN);
                                                                        if (!miss_ind[j]) { // never seen this indiv as missing before
                                                                            miss_ind[j] = true;
                                                                            n_miss_ind++;
                                                                        }
                                                                }
                                                                j++;
							}
							corr_index ++;
						}

						// Is this a valid lm model?
						if(lm -> isValid()){
							// One more valid test performing
							nepi ++;

							// Count as a good result
							summary_good[e1]++;
							if(sA[e2]) summary_good[e2] ++;

							// Do we want to record this as part of the summary for the first set?
							if (pvalue <= xpar::xepi_alpha2) {
								// first variable will always be in A set
								summary_sig[e1]++;

								// but the sechond may also be in A set
								if(sA[e2]) summary_sig[e2]++;
							}
						}

						// Is this result the best score yet for marker in set A
						if (pvalue <= best_score[e1]){
							best_score[e1] = pvalue;
							best_partner[e1] = e2;
						}

						// The second marker might also be in set A
						if(sA[e2]){
							if(pvalue <= best_score[e2]){
								best_score[e2] = pvalue;
								best_partner[e2] = e1;
							}
						}
										
						// Is this worth displaying?
						if(pvalue <= xpar::xepi_alpha1){
							XEPI << setw(4) << P.locus[e1]->chr << " "
							<< setw(par::pp_maxsnp) << P.locus[e1]->name << " "
							<< setw(4) << P.locus[e2]->chr << " "
							<< setw(par::pp_maxsnp) << P.locus[e2]->name << " ";

							if (lm -> isValid()){
								XEPI << setw(12) << b[3] << " "
								<< setw(12) << t_statistic_inter << " "
								<< setw(18) << pvalue << " "
								<< "\n";
							}else {
								XEPI << setw(12) << "NA" << " "
								<< setw(12) << "NA" << " "
								<< setw(12) << "NA" << " "
								<< "\n";
							}
							XEPI.flush();
						}
					} 
					else {
						cout << "Multicollinearity for " << P.locus[e1]->chr << " " << P.locus[e1]->name << " " << P.locus[e2]->chr << " " << P.locus[e2]->name << " (skipping)" << endl;
					}
					delete lm;

                                }// Next pair of SNPs
			}
		}
	}
	XEPI.close();
			
	//////////////////////
	// Summary of results
        XEPI.open((f + ".summary").c_str(), ios::out);
        XEPI.clear();
        
        P.printLOG("Performed a total of "+int2str(nepi)+" valid SNP X SNP tests\n");
        P.printLOG("Writing xepistasis summary results to [ " + f + ".summary ]\n");
        
        XEPI.precision(4);
        XEPI << setw(4) << "CHR" << " "
             << setw(par::pp_maxsnp) << "SNP" << " "
             << setw(12) << "N_SIG" << " "
             << setw(12) << "N_TOT" << " "
             << setw(12) << "PROP" << " "
             << setw(12) << "BEST_pvalue" << " "
             << setw(12) << "BEST_CHR" << " "
             << setw(par::pp_maxsnp) << "BEST_SNP" << " "
             << "\n";

        int c=0;
        for (int e1=0;e1<P.nl_all;e1++){
                if(sA[e1]){
                        XEPI << setw(4) << P.locus[e1]->chr << " "
                             << setw(par::pp_maxsnp) << P.locus[e1]->name << " "
                             << setw(12) << summary_sig[e1] << " "
                             << setw(12) << summary_good[e1] << " "
                             << setw(12) << (double)summary_sig[e1] / (double) summary_good[e1] << " "
                             << setw(12) << best_score[e1] << " "
                             << setw(4) << P.locus[best_partner[e1]]->chr << " "
                             << setw(par::pp_maxsnp) << P.locus[best_partner[e1]]->name << " "
                             << "\n";
                }
        }
	XEPI.close();

	////////////////////////////////////////////////////
	// Output correlation matrix (only for quantitative)
	if(par::qt){	
		f += ".corr";
		XEPI.open(f.c_str(), ios::out);
		XEPI.clear();
		XEPI.precision(12);

		P.printLOG("Writing correlation matrix for each pair of interactions to [ " + f + " ]\n");
		if(corr_sparse.size() == 0 || corr_sparse[0].size() == 0) {
			return;
		}

		matrix_t corr_pruned;
		matrix_t corr_pruned_t;
                sizeMatrix(corr_pruned,corr_index,P.n - n_miss_ind); // # snps X # indivs
                sizeMatrix(corr_pruned_t,P.n - n_miss_ind,corr_index); // # indivs X # snps

                for (int i=0;i<corr_index; i++){
                        int l = 0;
                        for (int j=0;j<P.n;j++){
                                if (!miss_ind[j]) {
                                        corr_pruned[i][l] = corr_sparse[i][j];
                                        corr_pruned_t[l][i] = corr_sparse[i][j]; // transpose
                                        l++;
                                }
                        }
                }

		matrix_t corr_matrix;
                multMatrix(corr_pruned,corr_pruned_t,corr_matrix);
	    		
                for(int i=0;i<corr_index;i++){
                        corr_matrix[i][i] = 1;
                        for( int j=0; j<corr_index-1; j++) {
                                XEPI << setw(12) << corr_matrix[i][j] << "\t";
                        }
                        XEPI << setw(12) << corr_matrix[i][corr_index-1] << endl;
                        XEPI.flush();
                }
                if (corr_index > 1000) {
                    cout << "WARNING: correlation matrix is too large to load into R." << endl;
                    cout << "You will not be able to run gene_based_interaction.sh." << endl;
                    cout << "For a smaller correlation matrix, consider trimming your sample size via LD pruning using the following commands:" << endl;
                    cout << "xwas --bfile [FILENAME] --indep-pairwise 50 5 0.3" << endl;
                    cout << "xwas --bfile [FILENAME] --extract xwas.prune.in --out [FILENAME_PRUNED]" << endl;
                    cout << "See PLINK documentation for full details." << endl;
                }
	}
	XEPI.close();
}
