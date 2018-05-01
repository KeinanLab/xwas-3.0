
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

void xStratSexAssocTest(Plink& P, matrix_t& resultpv_male, matrix_t& resultpv_female, vector<double>& pval){
	string f = par::output_file_name;
	//temp file for gene based testing added by Liang
	string f_stouffers = f+"_stouffer_temp.txt";
	string f_fishers = f+"_fisher_temp.txt";
	if (par::bt){
		f += ".xstrat.logistic";
		P.printLOG("Writing sex-stratified logistic model association results to [ " + f + " ] \n");
	} 
	else{
		f += ".xstrat.linear";
		P.printLOG("Writing sex-stratified linear model association results to [ " + f + " ] \n");
	}
	ofstream ASC;
	ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " "
	<< setw(par::pp_maxsnp) << "SNP" << " "
	<< setw(10) << "BP" << " "
	<< setw(4) << "A1" << " "
	<< setw(4) << "A2" << " "
	<< setw(10) << "TEST" << " ";
	ofstream ASC1;
	ofstream ASC2;

	if (par::bt && !xpar::xreturn_beta) ASC << setw(10) << "OR_M" << " ";
	else ASC << setw(10) << "BETA_M" << " ";
	ASC << setw(12) << "P_M" << " ";

	if (par::display_ci || par::qt) {
	  ASC << setw(12) << "SE_M" << " ";
	  if (par::bt)
	    ASC << setw(12) << string("L"+dbl2str(par::ci_level*100)+"_M") << " "
		<< setw(12) << string("U"+dbl2str(par::ci_level*100)+"_M") << " ";
	}

	if (par::bt && !xpar::xreturn_beta) ASC << setw(10) << "OR_F" << " ";
	else ASC << setw(10) << "BETA_F" << " ";
	ASC << setw(12) << "P_F" << " ";

	if (par::display_ci || par::qt) {
	  ASC << setw(10) << "SE_F" << " ";
	  if (par::bt)
	    ASC << setw(12) << string("L"+dbl2str(par::ci_level*100)+"_F") << " "
		<< setw(12) << string("U"+dbl2str(par::ci_level*100)+"_F") << " ";
	}
	
	// If-cases for determining headers
	if (xpar::fishers && !xpar::stouffers){
		ASC << setw(12) << "Fisher_Chi_Squared" << " ";
		ASC << setw(12) << "P_comb_Fisher" << " " << "\n";
		//create ASC1 for gene based testing
		ASC1.open(f_fishers.c_str(), ios::out); //file with snp names [column1] and p values [column2]
	}
	if (!xpar::fishers && xpar::stouffers){
		ASC << setw(12) << "Stouffers_Z" << " ";
		ASC << setw(12) << "P_comb_Stouffer" << " " << "\n";
		//create ASC1 for gene based testing
		ASC1.open(f_stouffers.c_str(), ios::out); //file with snp names [column1] and p values [column2]
	}
	if (xpar::fishers && xpar::stouffers){
		ASC << setw(12) << "Fisher_Chi_Squared" << " " 
		<< setw(12) << "P_comb_Fisher" << " ";
		ASC << setw(12) << "Stouffers_Z" << " "
		<< setw(12) << "P_comb_Stouffer" << " " << "\n";
		//create ASC1 for gene based testing
		ASC1.open(f_fishers.c_str(), ios::out); //file with snp names [column1] and p values [column2]
		ASC2.open(f_stouffers.c_str(), ios::out); 
	}
	ASC.precision(4);

	// combine the p-values
	vector<Locus*>::iterator loc = P.locus.begin();
	matrix_t::iterator pv_f_it = resultpv_female.begin();
	matrix_t::iterator pv_m_it = resultpv_male.begin();

	// Lauren added
	vector<CSNP*>::iterator s = P.SNP.begin();

	//only if there is something for output. If the user specifies --xchr-model 0, recessive, etc., it won't work for chrX.
	if ((resultpv_male.size() > 0) && (resultpv_female.size() > 0)){ 
	        // iterates over each locus (SNP)
		while ( loc != P.locus.end() ){
		        // Giant if statement: only prints info for a SNP if it passes this
		        // IF SNP IS ON X
			if (par::chr_sex[(*loc)->chr]){
			        // Added by Lauren 2/22/18
			        // Count alleles for bt SE
			        int mA1 = 0, mA2 = 0;
				int mU1 = 0, mU2 = 0;
			        int fA1 = 0, fA2 = 0;
				int fU1 = 0, fU2 = 0;

				// Statistics for qt SE
				int nanal_M = 0, nanal_F = 0;
				double qt_g_covar_M=0, qt_g_covar_F=0;
				double qt_mean_M=0, qt_var_M=0;
				double g_mean_M=0, g_var_M=0;
				double qt_mean_F=0, qt_var_F=0;
				double g_mean_F=0, g_var_F=0;

				/////////////////////////////
				// Iterate over individuals
				/////////////////////////////
      
				vector<bool>::iterator i1 = (*s)->one.begin();
				vector<bool>::iterator i2 = (*s)->two.begin();
				vector<Individual*>::iterator gperson = P.sample.begin();
 
				// Only iterate if need to calculate SE
				if (par::display_ci || par::qt) {
        				while ( gperson != P.sample.end() ) {
        				  // Phenotype for this person (i.e. might be permuted)
        				  Individual * pperson = (*gperson)->pperson;
        
        				  // Is this individual missing? Next person
        				  if ( pperson->missing ) {
        				    gperson++;
        				    i1++;
        				    i2++;
        				    continue;
        				  }
        	  
        				  // SNP alleles	  
        				  bool s1 = *i1;
        				  bool s2 = *i2;

					  // Add to num_in_analysis if non-missing genotype
					  if ( ( ! ( s1 && (!s2) ) ) && par::qt ) {
					    if ( (*gperson)->sex ) {
					      nanal_M++;
					      qt_mean_M += pperson->phenotype;
					    } else {
					      nanal_F++;
					      qt_mean_F += pperson->phenotype;
					    }
					  }

        				  // Sex determines if we add 1 allele or 2
        				  if ( (*gperson)->sex ) { // male
        				    if (pperson->aff) {  // affected 
        				      if (!s1) {
        					if (!s2) {
						  mA1++;
						  g_mean_M++;
						}
        				      }
        				      else {
        					if ( s2 ) mA2++;
        				      }
        				    } else { // unaffected
        				      if (!s1) {
        					if (!s2) mU1++;
        				      }
        				      else {
        					if ( s2) mU2++;
        				      }
        				    }
        				  } else { // female 
        				    if (pperson->aff) { // affected 
        				      if (!s1) {
        					if (!s2) {
						  fA1 += 2;
						  g_mean_F += 2;
						}
        					else {
        					  fA1++; fA2++;
						  g_mean_F++;
						}
        				      } else {
        					if ( s2 ) fA2 += 2;
        				      }
        				    } else {  // unaffected
        				      if (!s1) {
        					if (!s2) fU1 += 2;
        					else
        					  { fU1++; fU2++; } 
        				      }
        				      else {
        					if ( s2 ) fU2 += 2;
        				      }
        				    }
        				  }
        				  // Next person
        				  gperson++;
        				  i1++;
        				  i2++;
        				} // end while loop

					//////////////////////////////////
					// Iterate over individuals again

					if (par::qt) { // only if calculating SE for qt

					  qt_mean_M /= (double)nanal_M;
					  qt_mean_F /= (double)nanal_F;
					  g_mean_M /= (double)nanal_M;
					  g_mean_F /= (double)nanal_F;

					  i1 = (*s)->one.begin();
					  i2 = (*s)->two.begin();
					  gperson = P.sample.begin();
	        
					  while ( gperson != P.sample.end() ) {
	  	  
					    // Phenotype for this person (i.e. might be permuted)
					    Individual * pperson = (*gperson)->pperson;
	  
					    // Is this individual missing? Next person
					    if ( pperson->missing ) {
					      gperson++;
					      i1++;
					      i2++;
					      continue;
					    }
	  	  
					    // SNP alleles
					    bool s1 = *i1;
					    bool s2 = *i2;

					    double g = 0;

					    // If non-missing genotype
					    if ( ! ( s1 && (!s2) ) ) {
					      if ( (*gperson)->sex ) { // male
						qt_var_M += (pperson->phenotype - qt_mean_M) * ( pperson->phenotype - qt_mean_M );
						if (!s1) g = 1;
						g_var_M += (g - g_mean_M) * (g - g_mean_M);
						qt_g_covar_M += ( pperson->phenotype - qt_mean_M ) * ( g - g_mean_M );
					      } else { // female
						qt_var_F += (pperson->phenotype - qt_mean_F) * ( pperson->phenotype - qt_mean_F );
						if (!s1) {
						  if (!s2) g = 2;
						  else g = 1;
						}
						g_var_F += (g - g_mean_F) * (g - g_mean_F);
						qt_g_covar_F += ( pperson->phenotype - qt_mean_F ) * ( g - g_mean_F );
					      }
					    }
	  
					    // Advance to the next person
					    gperson++;
					    i1++;
					    i2++;
					  } // end while loop

					  qt_var_M /= (double)nanal_M - 1;
					  g_var_M /= (double)nanal_M - 1;
					  qt_g_covar_M /= (double)nanal_M - 1;

					  qt_var_F /= (double)nanal_F - 1;
					  g_var_F /= (double)nanal_F - 1;
					  qt_g_covar_F /= (double)nanal_F - 1;

					} // end if par::qt

				} // end if par::display_ci
        
			        // CHR, SNP, BP, A1, TEST
				ASC << setw(4) << (*loc)->chr << " "
				<< setw(par::pp_maxsnp) << (*loc)->name << " "
				<< setw(10) << (*loc)->bp << " "
				<< setw(4)  << (*loc)->allele1 << " "
				<< setw(4)  << (*loc)->allele2 << " "
				<< setw(10) << "SexStrat" << " ";

				double beta_m = *(pv_m_it->begin());
				double pvm = *(pv_m_it->begin() + 1);
				double beta_f = *(pv_f_it->begin());
				double pvf = *(pv_f_it->begin() + 1);
				double zero = 0;

				// OR_M, P_M
				if (par::bt && !xpar::xreturn_beta) { // binary trait and not display beta
				  if (pvm > 0) {
				    ASC << setw(10) << exp(beta_m) << " " 
					<< setw(12) << pvm << " ";
				  } else {
				    ASC << setw(10) << "NA" << " " 
					<< setw(12) << "NA" << " ";
				  }
				// BETA_M, P_M, SE_M
				} else { // non-binary trait or display beta
				  if (pvm > 0) {
				    ASC << setw(10) << beta_m << " " 
					<< setw(12) << pvm << " ";
				    if (par::qt) {
				      double vbeta_M = ( qt_var_M/g_var_M - (qt_g_covar_M*qt_g_covar_M)/(g_var_M*g_var_M) ) / (nanal_M-2);
				      ASC << setw(10) << sqrt(vbeta_M) << " ";
				    }
				  } else {
				    ASC << setw(10) << "NA" << " " 
					<< setw(12) << "NA" << " ";
				    if (par::qt) ASC << setw(10) << "NA" << " ";
				  }
				}

				// SE_M, L_M, U_M
				if (par::bt && par::display_ci) {
				  if (pvm > 0) {
				    double SE_M = sqrt(1/(double)mA1 + 1/(double)mA2 + 1/(double)mU1 + 1/(double)mU2);
				    double OR_M = exp(beta_m);
				    double lOR_M = log(OR_M);
				    double OR_lower_M = exp( lOR_M - par::ci_zt * SE_M );
				    double OR_upper_M = exp( lOR_M + par::ci_zt * SE_M );

				    if (SE_M != SE_M || SE_M == 1/zero || SE_M == -1/zero ) {
				      ASC << setw(12) << "NA" << " "
					  << setw(12) << "NA" << " "
					  << setw(12) << "NA" << " ";
				    } else {
				      ASC << setw(12) << SE_M << " "
					  << setw(12) << OR_lower_M << " "
					  << setw(12) << OR_upper_M << " ";
				    }
				  } else {
				    ASC << setw(12) << "NA" << " "
					<< setw(12) << "NA" << " "
					<< setw(12) << "NA" << " ";
				  }
				}

				// OR_F, P_F
				if (par::bt && !xpar::xreturn_beta) { // binary trait and not display beta
				  if (pvf > 0) {
				    ASC << setw(10) << exp(beta_f) << " " 
					<< setw(12) << pvf << " ";
				  }
				  else {
				    ASC << setw(10) << "NA" << " " 
					<< setw(12) << "NA" << " ";
				  }
				// BETA_F, P_F, SE_F
				} else { // non-binary trait or display beta
				  if (pvf > 0) {
				    ASC << setw(10) << beta_f << " " 
					<< setw(12) << pvf << " ";
				    if (par::qt) {
				      double vbeta_F = ( qt_var_F/g_var_F - (qt_g_covar_F*qt_g_covar_F)/(g_var_F*g_var_F) ) / (nanal_F-2);
				      ASC << setw(10) << sqrt(vbeta_F) << " ";
				    }
				  } else {
				    ASC << setw(10) << "NA" << " " 
					<< setw(12) << "NA" << " ";
				    if (par::qt) ASC << setw(10) << "NA" << " ";
				  }
				}

				// SE_F, L_F, U_F
				if (par::bt && par::display_ci) {
				  if (pvf > 0) {
				    double SE_F = sqrt(1/(double)fA1 + 1/(double)fA2 + 1/(double)fU1 + 1/(double)fU2);
				    double OR_F = exp(beta_f);
				    double lOR_F = log(OR_F);
				    double OR_lower_F = exp( lOR_F - par::ci_zt * SE_F );
				    double OR_upper_F = exp( lOR_F + par::ci_zt * SE_F );

				    if (SE_F != SE_F || SE_F == 1/zero || SE_F == -1/zero ) {
				      ASC << setw(12) << "NA" << " "
					  << setw(12) << "NA" << " "
					  << setw(12) << "NA" << " ";
				    } else {
				      ASC << setw(12) << SE_F << " "
					  << setw(12) << OR_lower_F << " "
					  << setw(12) << OR_upper_F << " ";
				    }
				  } else {
				    ASC << setw(12) << "NA" << " "
					<< setw(12) << "NA" << " "
					<< setw(12) << "NA" << " ";
				  }
				}

				double chisq = ((pvf > 0) && (pvm > 0)) ? -2*(log(pvf)+log(pvm)) : -1.;

				// We consider weighting by individual, not chromosome here
				int nmales = 0, nfemales = 0;
				for (int i=0; i<P.n; i++) {
					if ( ! P.sample[i]->missing ){
						if ( P.sample[i]->sex )	nmales++;
						else nfemales++;
					}
				}

				// Weighting the females more heavily (if the flag is there)
				if (xpar::sex_weight){
					nfemales = 2 * nfemales;
				}

				// Sample size weighting method (maybe implement inverse variance method sometime?)
				double wm = sqrt(nmales);
				double wf = sqrt(nfemales);	

				// Direction of effect from the odds ratio
				double efm = beta_m > 0 ? 1 : -1;
				double eff = beta_f > 0 ? 1 : -1;

				//Individual test statistics
				double zm = fabs(ltqnorm((pvm / 2.0))) * efm;
				double zf = fabs(ltqnorm((pvf / 2.0))) * eff;

				double stouffer_z = ((wm * zm) + (wf * zf)) / sqrt(nmales + nfemales);
				double pv_stouffer = 2 * (1 - normdist(fabs(stouffer_z)));
				//Weighting by Inverse variance (TODO)

				// Fisher_Chi_Squared, P_comb_Fisher
				if (xpar::fishers && !xpar::stouffers){
					if (chisq > 0)	{
   					        double pv_fishers = chiprobP(chisq,4);
						ASC << setw(12) << chisq << " " << setw(12) << pv_fishers;
						ASC1 << setw(par::pp_maxsnp) << (*loc)->name << " " << setw(12) << pv_fishers <<"\n";
						pval.push_back(pv_fishers);
					}
					else	ASC << setw(12) << "NA" << " " << setw(12) << "NA";
				}
				// Here we have set the parameters on this 
				// Stouffers_Z, P_comb_Stouffer
				if (!xpar::fishers && xpar::stouffers){
					if ((pvf > 0) && (pvf > 0)) {
						ASC << setw(12) << stouffer_z << " " << setw(12) << pv_stouffer;
						ASC1<< setw(par::pp_maxsnp) << (*loc)->name << " "<<setw(12)<< pv_stouffer<<"\n";
						pval.push_back(pv_stouffer);
					}
			  		else ASC << setw(12) << "NA" << " " << setw(12) << "NA";
				}
				// Fisher_Chi_Squared, P_comb_Fisher, Stouffers_Z, P_comb_Stouffer
				if (xpar::fishers && xpar::stouffers){
					// Need to establish other field for stouffers z-score
					if (chisq > 0) {
						ASC << setw(12) << chisq << " " << setw(12) << chiprobP(chisq, 4);
						ASC1<< setw(par::pp_maxsnp) << (*loc)->name << " "<<setw(12)<< chiprobP(chisq,4)<<"\n";
					}
					else ASC << setw(12) << "NA" << " " << setw(12) << "NA";
					if ((pvf > 0) && (pvf > 0)){
						ASC << setw(12) << stouffer_z << " " << setw(12) << pv_stouffer;
						ASC2<< setw(par::pp_maxsnp) << (*loc)->name << " "<<setw(12)<< pv_stouffer<<"\n";
					}
					else ASC << setw(12) << "NA" << " " << setw(12) << "NA";
				}
				ASC << "\n";
				pv_f_it ++;
				pv_m_it ++;
			}
			loc ++;
			s ++;
		}
	}
} // end xStratSexAssocTest
