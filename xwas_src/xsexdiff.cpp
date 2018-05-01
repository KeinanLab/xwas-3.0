
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

int makeRanks(vector< beta >& b) { //b is already sorted by value
	vector< beta >::iterator it = b.begin(), tieSt = it, tmpIt = it;
	double tmpValue = it -> value;
	int nTie = 1, tmpRank = 1, sumRanks = 1;
	it ++;
	while(1) {
		if ((it == b.end()) || (it -> value != tmpValue)) {
			double avgRank = 1. * sumRanks / nTie;
			for(tmpIt = tieSt; tmpIt != it; tmpIt ++) {
				tmpIt -> rank = avgRank;
			}
			if (it == b.end()) {
				break;
			}
			nTie = 1;
			tmpRank ++;
			tmpValue = it -> value;
			sumRanks = tmpRank;
			tieSt = it;
		}
		else { //tie case
			nTie ++;
			tmpRank ++;
			sumRanks += tmpRank;
		}
		it ++;
	}
	return 0; 
}

double corrRank(const vector< beta >& x, const vector< beta >& y) {
	double sumX = 0., sumY = 0., sumX2 = 0., sumY2 = 0., sumXY = 0., n = x.size();
	vector< beta >::const_iterator itX = x.begin(), itY = y.begin();
	for(; itX != x.end(); itX ++, itY ++) {
		sumX += itX -> rank;
		sumY += itY -> rank;
		sumX2 += square(itX -> rank);
		sumY2 += square(itY -> rank);
		sumXY += itX -> rank * itY -> rank;
	}
	//modified by Liang (consider the edge case)
	if (n*sumXY - sumX * sumY==0){
		return 0;
	}
	else{
		return (n * sumXY - sumX * sumY) / sqrt(n * sumX2 - square(sumX)) 
			/ sqrt(n * sumY2 - square(sumY));
		}
}

double spearman(vector< beta >& betaM, vector< beta >& betaF) {
	sort(betaM.begin(), betaM.end(), byValue());
	sort(betaF.begin(), betaF.end(), byValue());
	makeRanks(betaM); makeRanks(betaF);
	sort(betaM.begin(), betaM.end(), byIndex());
	sort(betaF.begin(), betaF.end(), byIndex());
	return corrRank(betaM, betaF);
}

double sexDiffPVal(double bM, double bF, double nM, double nF, double pM, double pF, double corr) { //using # of individuals
	double SEM = bM / (-ltqnorm(0.5 * pM));
	double SEF = bF / (-ltqnorm(0.5 * pF));  //modified by Liang
	double SEM2 = square(SEM), SEF2 = square(SEF);
	double SEMF = sqrt(SEM2 + SEF2 - 2. * corr * SEM * SEF);
	double stat = (bM - bF) / SEMF;
	double df = square(SEM2 / nM + SEF2 / nF) / (square(SEM2 / nM) / (nM - 1.) + 
				square(SEF2 / nF) / (nF - 1.));
	return calc_tprob(stat, df);
}

//sex difference test
void xSexDiffTest(Plink& P, matrix_t& resultpv_male, matrix_t& resultpv_female){
	string f = par::output_file_name;
	if (par::bt){
		f += ".xdiff.logistic";
		P.printLOG("Writing logistic model sex difference test results to [ " + f + " ] \n");
	} 
	else{
		f += ".xdiff.linear";
		P.printLOG("Writing linear model sex difference test results to [ " + f + " ] \n");
	}
	ofstream ASC;
	ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " "
	<< setw(par::pp_maxsnp) << "SNP" << " "
	<< setw(10) << "BP" << " "
	<< setw(4) << "A1" << " "
	<< setw(10) << "TEST" << " ";

	// We will actually need to start the parsing from here

	if (par::bt && !xpar::xreturn_beta)	ASC << setw(10) << "OR_M" << " ";
	else ASC << setw(10) << "BETA_M" << " ";
	ASC << setw(12) << "P_M" << " ";
	if (par::bt && !xpar::xreturn_beta)	ASC << setw(10) << "OR_F" << " ";
	else ASC << setw(10) << "BETA_F" << " ";
	ASC << setw(12) << "P_F" << " ";
	ASC << setw(12) << "P_DIFF" << "\n";
	
	ASC.precision(4);

	// combine the p-values
	vector<Locus*>::iterator loc = P.locus.begin();
	matrix_t::iterator pv_f_it = resultpv_female.begin();
	matrix_t::iterator pv_m_it = resultpv_male.begin();

	//only if there is something for output. If the user specifies --xchr-model 0, recessive, etc., it won't work for chrX.
	if ((resultpv_male.size() > 0) && (resultpv_female.size() > 0)){ 
		int n = 0;
		vector< beta > betaF, betaM;
		while ( loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){
				if ((*(pv_m_it->begin() + 1) > 0) && (*(pv_f_it->begin() + 1) > 0)) {
					beta tBM, tBF;
					tBM.index = n; tBF.index = n;
					tBM.value = *(pv_m_it->begin()); tBF.value = *(pv_f_it->begin());
					betaM.push_back(tBM);
					betaF.push_back(tBF);
					n ++;
				} 
				pv_f_it ++;
				pv_m_it ++;
			}
			loc ++;
		}
		double r = spearman(betaM, betaF);
		loc = P.locus.begin();
		pv_f_it = resultpv_female.begin();
		pv_m_it = resultpv_male.begin();
		while ( loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){
				ASC << setw(4)  << (*loc)->chr << " "
				<< setw(par::pp_maxsnp) << (*loc)->name << " "
				<< setw(10) << (*loc)->bp << " "
				<< setw(4)  << (*loc)->allele1 << " "
				<< setw(10) << "SexDiff" << " ";
				if (par::bt && !xpar::xreturn_beta){
					if (*(pv_m_it->begin() + 1) > 0)	
						ASC << setw(10) << exp(*(pv_m_it->begin())) << " " 
						<< setw(12) << *(pv_m_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
					if (*(pv_f_it->begin() + 1) > 0)
						ASC << setw(10) << exp(*(pv_f_it->begin())) << " " 
						<< setw(12) << *(pv_f_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
				}
				else{
					if (*(pv_m_it->begin() + 1) > 0)	
						ASC << setw(10) << *(pv_m_it->begin()) << " " 
						<< setw(12) << *(pv_m_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
					if (*(pv_f_it->begin() + 1) > 0)
						ASC << setw(10) << *(pv_f_it->begin()) << " " 
						<< setw(12) << *(pv_f_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
				}
				
				// We consider weighting by individual, not chromosome here
				int nmales = 0, nfemales = 0;
				for (int i=0; i<P.n; i++) {
					if ( ! P.sample[i]->missing ){
						if ( P.sample[i]->sex )	nmales++;
						else nfemales++;
					}
				}

				//should we initiate some check in here? to make sure they are below 1?
				if ((*(pv_m_it->begin() + 1) > 0) && (*(pv_f_it->begin() + 1) > 0) && (*(pv_m_it->begin() + 1) < 1) && (*(pv_f_it->begin() + 1) < 1)) {
					double bM = *(pv_m_it->begin()), bF = *(pv_f_it->begin());
					double pM = *(pv_m_it->begin() + 1), pF = *(pv_f_it->begin() + 1);
					ASC << setw(12) << sexDiffPVal(bM, bF, nmales, nfemales, pM, pF, r) << " ";
				}
				else {
					ASC << setw(12) << "NA" << " ";
				}
				
				ASC << "\n";
				pv_f_it ++;
				pv_m_it ++;
			}
			loc ++;
		}
	}
}
