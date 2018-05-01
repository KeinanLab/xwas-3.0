
///////////////////////////////////////////////////////////////////////////
//       Edited by Liang Zhang (SNP level clayton's test)                //
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

void xClayton (Plink& P, vector<double>& pval_results){
	string f = par:: output_file_name;
	//temp file for gene_based testing
	string f1 = f +"_temp.xclayton";
	if (par::qt){
		f +=".xclayton";
		P.printLOG("Writing Clayton's test results to [ " + f + " ] \n");
		ofstream ASC;
		ASC.open(f.c_str(), ios::out);
		ASC << setw(4) << "CHR" << " "
		<< setw (par::pp_maxsnp) << "SNP" << " "
		<< setw(10) << "BP" << " "
		<< setw(4) << "A1" << " "
		<< setw(10) << "CHISQ-STAT" << " "
		<< setw(10) << "DF" << " "
		<< setw(10) << "PVAL\n";
		//create ASC1 for gene_based testing
		ofstream ASC1;
		ASC1.open(f1.c_str(), ios::out); //file with snp names [column1] and p-values [column2]

		int nind = P.sample.size(); //total number of individuals

		vector<Locus*>::iterator loc = P.locus.begin();
		for (int l=0; l < P.locus.size();l++){
		        // CHR, SNP, BP, A1
			ASC << setw(4) << (*loc)->chr << " "
			<< setw(par::pp_maxsnp) << (*loc)->name << " "
			<< setw(10) << (*loc)->bp << " "
			<< setw(4)  << (*loc)->allele1 << " ";

			//need to recalculate all these in order to consider missingness for an individual
			vector_t Y, YF;
			int nfemales=0;
			int nmales=0;
			// Iterate over each individual
			for (int i=0;i<nind; i++){
				Individual * person = P.sample[i];
				if ((!person->missing) && (!MISSING(i,l))) {
					if (person->sex){
						nmales++;
					}
					else{
						nfemales++;
						YF.push_back(person->phenotype);
					}
					Y.push_back(person->phenotype);
				}
			}
			double meanY=calcMean(Y);
			double meanYF=calcMean(YF);

			double U_D=0.;
			double U_A=0.;
			double D_f=0.;

			vector_t Xa(nind), Xd(nind);
			for (int i=0; i< nind; i++){
				Individual * person = P.sample[i];
				if(!(person->missing) && !(MISSING(i,l))){  //not missing
					int xa = 1;
					int xd = 1;
					bool s1 = P.SNP[l]->one[i];
					bool s2 = P.SNP[l]->two[i];
					if (s1 && s2){ 		//both are major (homozogous) 11 = hom(22)
						xa = 2;
						xd = 0;
					}
					else if (!s1 && !s2){	//both are minor (homozogous) 00 = hom(11)
						xa = 0;
						xd = 0;
					}
					U_A += (person->phenotype - meanY) * xa;
					if (!(person->sex)){  //only consider D for females
						U_D += (person->phenotype - meanYF) * xd;
						D_f += xd;
					}
					Xa[i] = xa;
					Xd[i] = xd;
				}
			}


			//since we assume allele freq are equal between males and females, 
			//meanA can be calculated form entire sample
			double meanXa=calcMean(Xa);
			double meanXd=double(D_f)/nfemales; //only consider females' D
			double sqAdiff=0.;
			double adint=0.;
			double sqDdiff=0.;
			double sqYF=0.;
			double sqYM=0.;

			for (int i=0; i< nind; i++){
				Individual * person = P.sample[i];
				if (!(person->missing) && !MISSING(i,l)){
					if (!(person->sex)){  //not male
						sqAdiff += square(Xa[i]-meanXa);
						sqDdiff += square(Xd[i]-meanXd);
						adint += (Xa[i]-meanXa)*(Xd[i]-meanXd);
						sqYF += square(person->phenotype - meanY);
					}
					else{
						sqYM += square(person->phenotype - meanY);	
					}
				}

			}

			matrix_t VF, VM, V;
			sizeMatrix(VF, 2, 2);
			sizeMatrix(VM, 2, 2);
			sizeMatrix(V, 2, 2);

			//scaled VF, already taken into consideration of sqYF
			VF[0][0] = sqAdiff * sqYF/(nfemales-1);
			VF[0][1] = adint * sqYF/(nfemales-1);
			VF[1][0] = adint * sqYF/(nfemales-1);
			VF[1][1] = sqDdiff * sqYF/(nfemales-1);
			//scaled VM, already taken into consideration of sqYM
			double p=meanXa/2;
			VM[0][0] = 4 * p * (1 - p) * sqYM;
			VM[0][1] = VM[1][0] = VM[1][1] = 0.;
			//V = VF+VM
			V[0][0]= VF[0][0]+VM[0][0];
			V[0][1]= VF[0][1]+VM[0][1];
			V[1][0]= VF[1][0]+VM[1][0];
			V[1][1]= VF[1][1]+VM[1][1];

			matrix_t V_inv;
			V_inv = inverse(V);
			matrix_t Ut;
			matrix_t U;
			sizeMatrix(Ut, 1, 2);
			sizeMatrix(U, 2, 1);
			U[0][0]= U_A;
			U[1][0]= U_D;
			Ut[0][0] = U_A;
			Ut[0][1] = U_D;
			matrix_t tmp;
			matrix_t chisq_matrix;
			multMatrix(Ut, V_inv, tmp);
			multMatrix(tmp, U, chisq_matrix);
			double chisq_r=chisq_matrix[0][0]; 
			double pval = chiprobP(chisq_r, 2);
			pval_results.push_back(pval);

			// CHISQ-STAT
			ASC << setw(10);
			if ( chisq_r == 0 ) ASC << "INF";
			else if( chisq_r < 0 ) ASC << "NA";
			else if( !realnum(chisq_r) ) ASC << "NA";
			else ASC << chisq_r; 
			ASC << " ";
			// DF
			ASC << setw(10) << 2 << " ";
			// PVAL
			ASC << setw(10);
			if ( pval == 0 ) ASC << "INF";
			else if( pval < 0 ) ASC << "NA";
			else ASC << pval; 
			ASC << "\n";

			//output temp file for gene based test
			ASC1 << setw(par::pp_maxsnp) << (*loc)->name <<" "<< setw(10) << pval << " "; 
			ASC1 <<"\n";
			++loc;
		}
	}else{
		error("--clayton-x for quantitative traits only");
	}
}
