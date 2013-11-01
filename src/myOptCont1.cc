#include"habitat4.h"

myOptContr::myOptContr(){};
myOptContr::~myOptContr(){};

void myOptContr::setVals( const SEXP &Rmaxit, const SEXP &RmaxitInner, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv, const SEXP &RGSiters)
{

	maxitQN = *(INTEGER( Rmaxit));
	maxitGS = *(INTEGER(RGSiters));
	maxitInner = *(INTEGER(RmaxitInner));
	traceQN = *(INTEGER(Rtrace));
	traceGS = 0;
	nReport = *(INTEGER(RnReport));
	abstol = *(REAL(Rabstol));
	reltol = *(REAL(Rreltol));
	reltolInner = reltol;
	conv = INTEGER( Rconv);

	denomEps = 1e-5;
}
