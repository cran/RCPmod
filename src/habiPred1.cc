#include"habitat4.h"
#include"habiPred1.h"

extern "C" { SEXP Habitat_predict_C( SEXP RX, SEXP Ry,
																	 SEXP Ralpha, SEXP Rtau, SEXP Rbeta,
																	 SEXP RalphaBoot, SEXP RtauBoot, SEXP RbetaBoot,
																	 SEXP RS, SEXP RG, SEXP RnObs, SEXP Rp, SEXP Rnboot,
																	 SEXP RptPreds, SEXP RbootPreds, SEXP Rtype,
																	 SEXP Rconc, SEXP Rsd)
{
	allClasses all;
	int nboot = *(INTEGER(Rnboot));
	double *bootalpha, *boottau, *bootbeta, *bootParms;
	double *ptPreds, *bootPreds;
	int type = *(INTEGER( Rtype));
	int bootCount;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals( Ry, RX, RS, RG, Rp, RnObs);	//read in the data
	all.parms.setVals( all.data, Ralpha, Rbeta, Rtau, Rconc, Rsd);	//read in the parameters

	vector<double> logPis(all.data.nG, all.data.NAnum);//, pis( dat.nG, dat.NAnum);
	vector< vector<double> > allPis( all.data.nObs, vector<double> (all.data.nG, all.data.NAnum));	//(nObs x nG) matrix --2D
	vector<double> allMus( all.data.nG*all.data.nS, all.data.NAnum);	//(nG x nS) matrix

	ptPreds = REAL( RptPreds);
	for( int i=0; i<all.data.nObs; i++)
		calcLogPis( logPis, allPis.at(i), all.data, all.parms, i);
	if( type == 1){
		for( int g=0; g<all.data.nG; g++)
			for( int i=0; i<all.data.nObs; i++)
				ptPreds[MATREF(i,g,all.data.nObs)] = allPis.at(i).at(g);
	}
	if( type == 0){
		calcMuFits( allMus, all.data, all.parms);
		bootCount = 0;
		calcMargFits( ptPreds, bootCount, allMus, allPis, all.data);
	}

	//setting up the bootstrap values for alpha, tau and beta
	bootalpha = REAL( RalphaBoot);
	boottau = REAL( RtauBoot);
	bootbeta = REAL( RbetaBoot);
	bootParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	int kount;

	bootPreds = REAL( RbootPreds);
	for( int b=0; b<nboot; b++){
		kount = 0;
		for( int i=0; i<all.parms.nalpha; i++){
			bootParms[kount] = bootalpha[MATREF(b,i,nboot)];
			kount++;
		}
		for( int i=0; i<all.parms.ntau; i++){
			bootParms[kount] = boottau[MATREF(b,i,nboot)];
			kount++;
		}
		for( int i=0; i<all.parms.nbeta; i++){
			bootParms[kount] = bootbeta[MATREF(b,i,nboot)];
			kount++;
		}
		all.parms.update( bootParms, all.data);

		//calculate pis and store them
		int place;
		for( int i=0; i<all.data.nObs; i++)
			calcLogPis( logPis, allPis.at(i), all.data, all.parms, i);
		if( type == 1){
			for( int g=0; g<all.data.nG; g++)
				for( int i=0; i<all.data.nObs; i++){
					place = MATREF( MATREF(i,g,all.data.nObs),b,(all.data.nObs*all.data.nG));
					//place = MATREF( MATREF(i,all.data.nObs*all.data.nS+g,all.data.nObs),b,(all.data.nObs*(all.data.nS+all.data.nG)));
					bootPreds[place] = allPis.at(i).at(g);
				}
		}
		if( type == 0){
			calcMuFits( allMus, all.data, all.parms);
			calcMargFits( bootPreds, b, allMus, allPis, all.data);
		}

	}

	SEXP Rres;	//R object to return -- it is meaningless!
	Rres = allocVector(REALSXP,1);
	double *R_pointer = REAL(Rres);
	R_pointer[0] = -9999;

	return( Rres);


}

}

void calcMargFits( double *ptPreds, int bootCount, const vector<double> &allMus, const vector< vector<double> > &allPis, const myData &dat)
{
	for( int i=0; i<dat.nObs; i++)
		for( int s=0; s<dat.nS; s++){
			ptPreds[MATREF(MATREF(i,s,dat.nObs),bootCount,(dat.nObs*dat.nS))] = 0.0;
			for( int g=0; g<dat.nG; g++)
				ptPreds[MATREF(MATREF(i,s,dat.nObs),bootCount,(dat.nObs*dat.nS))] += allMus.at(MATREF(g,s,dat.nG))*allPis.at(i).at(g);
		}

}
