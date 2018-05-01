
///////////////////////////////////////////////////////////////////////////
//       Genomic control functionality: added by Lauren, 3/1/18          //
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

void xMultComp(Plink& P, vector<double>& pval, string f) {
  if (pval.size() == 0) return;

  // Convert pvals to chi-sq with df = 1
  vector<double> chi;
  vector<double>::iterator pval_it = pval.begin();
  while (pval_it != pval.end()) {
    double chi_sq = inverse_chiprob(*pval_it,1);
    chi.push_back(chi_sq);
    pval_it++;
  }

  double lambda = median(chi) / 0.456;
  if (lambda < 1) lambda = 1.00;

  P.printLOG("Genomic inflation factor (based on median chi-squared) is " + dbl2str(lambda)+"\n");
  P.printLOG("Correcting for " + int2str(chi.size()) + " tests\n");

  ofstream ASC;
  f += ".adjusted";
  ASC.open(f.c_str(),ios::out);
  ASC.precision(4);
  
  P.printLOG("Writing multiple-test corrected significance values to [ " + f + " ] \n");

  ASC << setw(4) << "CHR" << " "
     << setw(par::pp_maxsnp) << "SNP" << " "
     << setw(10) << "UNADJ" << " "
     << setw(10) << "GC" << " \n";

  vector<Locus*>::iterator loc = P.locus.begin();
  vector<double>::iterator chi_it = chi.begin();

  /////////////////////////////
  // Iterate over SNPs
  /////////////////////////////

  while ( loc != P.locus.end() ) {
    if (par::chr_sex[(*loc)->chr]) { // if on X

      // CHR, SNP
      ASC << setw(4) << (*loc)->chr << " "
	  << setw(par::pp_maxsnp) << (*loc)->name << " ";

      // Unadjusted
      double p = chiprobP(*chi_it,1);
      ASC << setw(10);
      if ( p == 0 ) ASC << "INF";
      else if( p < 0 ) ASC << "NA";
      else ASC << p; 
      ASC << " ";

      // Genomic control
      double gc = chiprobP(*chi_it / lambda, 1);
      ASC << setw(10);
      if ( gc == 0 ) ASC << "INF";
      else if( gc < 0 ) ASC << "NA";
      else ASC << gc;
      ASC << "\n";

      chi_it++;
    }
    loc++;
  }
  ASC.close();  
}
