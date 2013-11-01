#include"habitat4.h"

myData::myData(){};
myData::~myData(){};

void myData::setVals( SEXP &Ry, SEXP &RX, const SEXP &RS, const SEXP &RG, const SEXP &Rp, const SEXP &RnObs){

//	double *tmpD;

	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	np = *(INTEGER( Rp));
	nObs = *(INTEGER( RnObs));
	NAnum = -999999;

/*	tmpD = REAL( Ry); y.assign( tmpD, tmpD + LENGTH( Ry));
	tmpD = REAL( RX); X.assign( tmpD, tmpD + LENGTH( RX));*/
	y = REAL( Ry);
	X = REAL( RX);
}
