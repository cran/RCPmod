#include"habitat4.h"

myDerivs::myDerivs(){};
myDerivs::~myDerivs(){};

void myDerivs::setVals( const myData &dat, SEXP &RderivsAlpha, SEXP &RderivsTau, SEXP &RderivsBeta, SEXP &Rscores)
{
	Alpha = REAL( RderivsAlpha);
	Tau = REAL( RderivsTau);
	Beta = REAL( RderivsBeta);
	Scores = REAL( Rscores);
}

void myDerivs::zeroDerivs( const myData &dat)
{
	for( int i=0; i<dat.nS; i++)
		Alpha[i] = 0.0;
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		Tau[i] = 0.0;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++)
		Beta[i] = 0.0;
}

void myDerivs::updateDerivs( const myData &dat, const vector<double> &alphaDerivsI, const vector<double> &tauDerivsI, const vector<double> &betaDerivsI)
{
	for( int s=0; s<dat.nS; s++)
		Alpha[s] += alphaDerivsI.at(s);
	for( int g=0; g<(dat.nG-1); g++)
		for( int s=0; s<dat.nS; s++)
			Tau[MATREF(g,s,(dat.nG-1))] += tauDerivsI.at(MATREF(g,s,(dat.nG-1)));
	for( int g=0; g<(dat.nG-1); g++)
		for( int p=0; p<dat.np; p++)
			Beta[MATREF(g,p,(dat.nG-1))] += betaDerivsI.at(MATREF(g,p,(dat.nG-1)));
}

void myDerivs::updateScores( const myData &dat, const vector<double> &alphaDerivsI, const vector<double> &tauDerivsI, const vector<double> &betaDerivsI, int i)
{
	for( int s=0; s<dat.nS; s++)
		Scores[MATREF(i,s,dat.nObs)] = alphaDerivsI.at(s);
	for( int g=0; g<(dat.nG-1); g++)
		for( int s=0; s<dat.nS; s++)
			Scores[MATREF(i,(dat.nS+MATREF(g,s,(dat.nG-1))),dat.nObs)] = tauDerivsI.at(MATREF(g,s,(dat.nG-1)));
	for( int g=0; g<(dat.nG-1); g++)
		for( int p=0; p<dat.np; p++)
			Scores[MATREF(i,(dat.nS+(dat.nG-1)*dat.nS+MATREF(g,p,(dat.nG-1))),dat.nObs)] = betaDerivsI.at(MATREF(g,p,(dat.nG-1)));
}

void myDerivs::update( double *grArr, const myData &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		Tau[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		Beta[i] = grArr[kount];
		kount++;
	}
}

void myDerivs::getArray( double *grArr, const myData &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		grArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		grArr[kount] = Tau[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		grArr[kount] = Beta[i];
		kount++;
	}

}
