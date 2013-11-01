#ifndef habiPred1_hh
#define habiPred1_hh

#include<R.h>
#include<Rmath.h>
#include<Rinternals.h>
#include<R_ext/Applic.h>
#include<vector>

#undef length
#include <iostream>

using namespace std;
using std::vector;         // use vector as abbreviation for std::vector
using std::cout;

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////////
////////////////////////////////////////////////////////
extern "C" SEXP Habitat_predict_C( SEXP RX, SEXP Ry,
															SEXP Ralpha, SEXP Rtau, SEXP Rbeta,
															SEXP RalphaBoot, SEXP RtauBoot, SEXP RbetaBoot,
															SEXP RS, SEXP RG, SEXP RnObs, SEXP Rp, SEXP nboot,
															SEXP RptPreds, SEXP bootPreds, SEXP Rtype,
															SEXP Rconc, SEXP Rsd);
void calcMargFits( double *ptPreds, int bootCount, const vector<double> &allMus, const vector< vector<double> > &allPis, const myData &dat);


#endif
