#include"habitat4.h"

myParms::myParms(){};
myParms::~myParms(){};

void myParms::getAllTaus( vector<double> &newTau, const myData &dat) const
{
	double su;

	newTau.assign(dat.nG*dat.nS, dat.NAnum);
	//calculate sum-to-zero taus
	for( int s=0; s<dat.nS; s++){
		su = 0.0;
		for( int g=0; g<(dat.nG-1); g++){
			newTau.at( MATREF(g,s,dat.nG)) = Tau[MATREF(g,s,(dat.nG-1))];
			su += Tau[MATREF(g,s,(dat.nG-1))];
		}
		newTau.at( MATREF((dat.nG-1),s,dat.nG)) = -su;
	}

}

void myParms::setVals( const myData &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rtau, SEXP &Rconc, SEXP &Rsd)
{
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Tau = REAL( Rtau);
	Beta = REAL( Rbeta);
	conc = *(REAL( Rconc));
	sd = *(REAL( Rsd));

	nalpha = dat.nS;
	ntau = (dat.nG-1)*dat.nS;
	nbeta = (dat.nG-1)*dat.np;
	nTot = nalpha + ntau + nbeta;

}

void myParms::getArray(double *parArr, const myData &dat) const
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		parArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		parArr[kount] = Tau[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		parArr[kount] = Beta[i];
		kount++;
	}
}

void myParms::update( double *parArr, const myData &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		Tau[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		Beta[i] = parArr[kount];
		kount++;
	}
}
