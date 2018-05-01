

///////////////////////////////////////////////////////////////////////////
//       Adapted from plink source codes by Purcell et al. 2007          //
//       By Feng Gao and Diana Chang                                     //
///////////////////////////////////////////////////////////////////////////


#ifndef __XOPTIONS_H__
#define __XOPTIONS_H__

#include "plink.h"

class xpar 
{
 public:
  static bool xwas; //use xwas?
  static bool strat_sex; //stratify sex?
  static bool multi_xchr_model; // code males 0/1 0/2 concurrently
  static bool fishers; // stratify sex test with fishers test
  static bool stouffers; // stratify sex with stouffers test
  static bool sex_weight; // if this is on, we will weight females as twice
  static bool af_sex; //output sex freq report?
  static bool am_sex_test; //missingness difference between sexes?
  static double am_sex_test_limit; //p-value for missingness between sexes
  static bool af_sex_test; //allele frequency between sexes?
  static double af_sex_test_limit; //p-value for af_sex_test
  static bool xreturn_beta; //return beta for logistic regression
  static bool sex_diff; //sex difference test?
  static bool var_het; //variance heterogenity test?
  static bool xepi; //calculate SNPxSNP interaction
  static double xepi_alpha1; //p-value threshold for retain interaction pairs
  static double xepi_alpha2; //p-value for retain interaction for each SNP in summary file
  static bool Set_by_Set; // setting option in xepi
  static bool var_het_weight; //weighted variance heterogeneity test
  static bool var_het_comb; //combined variance heterogeneity test
  static bool clayton_x; //clayton's test on sex chromosome?
};

struct beta { 
  double value;
  size_t index;
  double rank;
};

struct byValue { 
    bool operator()(beta const& b1, beta const& b2) { 
        return b1.value < b2.value;
    }
};

struct byIndex { 
    bool operator()(beta const& b1, beta const& b2) { 
        return b1.index < b2.index;
    }
};

void xSetOptions(CArgs &);
void xfilterSNPs(Plink &);
void xStratSexAssocTest(Plink &, matrix_t &, matrix_t &, vector<double> &);
void xSexDiffTest(Plink &, matrix_t &, matrix_t &);
void xSexRelatedTests(Plink &);
void xMultComp(Plink& P, vector<double>&, string);
void xClayton(Plink &, vector<double>&);
void xVarianceHeterogeneity(Plink& P, vector<double>&);
void xAssocFunctions(Plink &);
void xParseArgs(CArgs &);
void calcXepistasis(Plink &);

#endif
