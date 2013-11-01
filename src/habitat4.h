#ifndef habitat1_hh
#define habitat1_hh

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
/////////////	Inline Function Definitions	//////////////
////////////////////////////////////////////////////////

#define MATREF( i, j, nrow) i + j*nrow	//indices start at 0

////////////////////////////////////////////////////////
////////////	Class Definitions	////////////////////////
////////////////////////////////////////////////////////

class myData
{
	public:
		myData();
		~myData();
		void setVals( SEXP &Ry, SEXP &RX, const SEXP &RS, const SEXP &RG, const SEXP &Rp, const SEXP &RnObs);

		int np, 		//the number of parameters in each of the (G-1) habitat lps
				nG,			//the number of habitats
				nS, 		//the number of species
				nObs,		//the number of observations
				NAnum;	//a common number to insert for NAs

		double 	*X, //the design matrix in vector form (nObs x np)
						*y;//the outcome matrix, in vector form (nObs x nS)
};

class myParms
{
	public:
		myParms();
		~myParms();
		void setVals( const myData &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rtau, SEXP &Rconc, SEXP &Rsd);
		void getArray(double *parArr, const myData &dat) const;
		void update( double *parArr, const myData &dat);
		void getAllTaus( vector<double> &newTau, const myData &dat) const;

//		void updatePars( double *newPar);

		double 	*Alpha, //the species' prevalences (Gx1)==(1xG)
						*Tau, 	//the species' habitat adjustments ((G-1)xS) -- only free params!
						*Beta,	//the habitats' free covariate params ((G-1)xp)
						conc,		//the concentration for the dirichlet penalty for pi
						sd;			//the sd for the normal penalty for tau
		int nalpha, ntau, nbeta, nTot;
//		int *myMask;
};

class myDerivs
{
	public:
		myDerivs();
		~myDerivs();
		void setVals( const myData &dat, SEXP &RderivsAlpha, SEXP &RderivsTau, SEXP &RderivsBeta, SEXP &Rscores);
		void zeroDerivs( const myData &dat);
		void updateDerivs( const myData &dat, const vector<double> &alphaDerivsI, const vector<double> &tauDerivsI, const vector<double> &betaDerivsI);
		void updateScores( const myData &dat, const vector<double> &alphaDerivsI, const vector<double> &tauDerivsI, const vector<double> &betaDerivsI, int i);
		void getArray( double *grArr, const myData &dat);
		void update( double *grArr, const myData &dat);

		double 	*Alpha, //the derivatives of logl w.r.t. alpha
						*Tau, 	//the derivatives of logl w.r.t. tau
						*Beta,	//the derivatives of logl w.r.t. beta
						*Scores;	//the score contribution for each site
};

class myOptContr
{
	public:
		myOptContr();
		~myOptContr();
		void setVals( const SEXP &Rmaxit, const SEXP &RmaxitInner, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv, const SEXP &RGSiters);

		int maxitQN, maxitGS, maxitInner, traceQN, traceGS, nReport, fnKount, grKount, ifail, *conv;
		double abstol, reltol, denomEps, reltolInner;
};

class allClasses
{
	public:
	allClasses();
	~allClasses();

	myData data;
	myParms parms;
	myDerivs derivs;
	myOptContr contr;
};

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////////
////////////////////////////////////////////////////////
extern "C" SEXP HABITAT_C( SEXP Ry, SEXP RX,
														 SEXP RS, SEXP RG, SEXP Rp, SEXP RnObs,
														 SEXP Ralpha, SEXP Rtau, SEXP Rbeta, SEXP Rconc, SEXP Rsd,
														 SEXP RderivsAlpha, SEXP RderivsTau, SEXP RderivsBeta, SEXP Rscores,
														 SEXP Rpis, SEXP Rmus, SEXP RlogDens, SEXP Rlogli,
														 SEXP Rmaxit, SEXP RmaxitInner, SEXP Rtrace, SEXP RnReport, SEXP Rabstol,
														 SEXP Rreltol, SEXP Rconv, SEXP RGSiters,
														 SEXP RloglOnly);

double invLogit( double x);
double mixLogl( const myData &dat, const myParms &parms, vector< vector<double> > &allPis, vector<double> &allMus, vector< vector<double> > &allLogDens, vector<double> &allLogls);
double calcTauPen( const myData &dat, const myParms &parms);
double calcPiPen( const vector<double> &logPis, const myData &dat, const myParms &parms);
void calcLogPis( vector<double> &logPis, vector<double> &pis, const myData &dat, const myParms &parms, 		int i);
void calcLogCondDens( vector<double> &condDens, const vector<double> &fits, const myData &dat, const myParms &parms, int i);
void calcMuFits( vector<double> &fits, const myData &dat, const myParms &parms);
double logBernoulli( double y, double mu);
double calcMixSum( const vector<double> &logPis, const vector<double> &condDens, double &wi, 							vector<double> &wij, int &maxg);
void loglDerivs( const myData &dat, const myParms &parms, myDerivs &derivs);
void calcDerivMu( vector<double> &muDerivs, const vector<double> &fits, const myData &dat, const 					myParms &parms, const double wi, const vector<double> &wij, const int m, const int i);
void calcDerivEtaMu( vector<double> &etaDerivsI, const myData &dat, const vector<double> &muDerivsI, const vector<double> &fits);
double logBernDer( double y, double mu);
void calcAlphaDeriv( vector<double> &alphaDerivsI, const vector<double> &etaDerivs, const myData &dat);
void calcTauDeriv( vector<double> &tauDerivsI, const vector<double> &etaDerivs, const myData &dat, const myParms &parms);
void calcTauPenDeriv( vector<double> &tauDerivsI, const myData &dat, const myParms &parms);
void calcPiDeriv( vector<double> &piDerivs, const myData &dat, const myParms &parms, const vector<double> &pis, const double 		wi, const vector<double> &wij, int m);
void calcBetaDeriv( vector<double> &betaDerivsI, const vector<double> &piDerivsI, const vector<double> 		&pis, const myData &dat, int i);
double GSoptimise( allClasses &all);
double ALLoptimise( allClasses &all);
bool converged( double *oldP, double *newP, const myOptContr &contr, int nTot);
double optimise_function(int n, double *par, void *ex);
void gradient_function(int n, double *par, double *gr, void *ex);

#endif
