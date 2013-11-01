#include"habitat4.h"

extern "C" { SEXP HABITAT_C( SEXP Ry, SEXP RX,
														 SEXP RS, SEXP RG, SEXP Rp, SEXP RnObs,
														 SEXP Ralpha, SEXP Rtau, SEXP Rbeta, SEXP Rconc, SEXP Rsd,
														 SEXP RderivsAlpha, SEXP RderivsTau, SEXP RderivsBeta, SEXP Rscores,
														 SEXP Rpis, SEXP Rmus, SEXP RlogDens, SEXP Rlogli,
														 SEXP Rmaxit, SEXP RmaxitInner, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP RGSiters,
														 SEXP RloglOnly)
{
	allClasses all;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals( Ry, RX, RS, RG, Rp, RnObs);	//read in the data
	all.parms.setVals( all.data, Ralpha, Rbeta, Rtau, Rconc, Rsd);	//read in the parameters
	all.derivs.setVals( all.data, RderivsAlpha, RderivsTau, RderivsBeta, Rscores);
	all.contr.setVals( Rmaxit, RmaxitInner, Rtrace, RnReport, Rabstol, Rreltol, Rconv, RGSiters);

	vector< vector<double> > allPis( all.data.nObs, vector<double> (all.data.nG, all.data.NAnum));	//(nObs x nG) matrix --2D
	vector<double> allMus( all.data.nG*all.data.nS, all.data.NAnum);	//(nG x nS) matrix
	vector< vector<double> > allLogDens( all.data.nObs, vector<double> (all.data.nG, all.data.NAnum));
	vector<double> allLogls( all.data.nObs, all.data.NAnum);

	double logl, logl1;

	//doing the optimising
	logl = mixLogl( all.data, all.parms, allPis, allMus, allLogDens, allLogls);
	if( *REAL(RloglOnly) != 1){
//		cout.precision(5);
//		cout << "Initial logl: " << fixed << -logl << "\n";
		Rprintf( "Initial logl: %10.5f \n", -logl);
		if( all.contr.maxitGS > 0)
			logl1 = GSoptimise( all);
		logl = ALLoptimise( all);
	}
	else
		loglDerivs( all.data, all.parms, all.derivs);

	//re-running to get pis and mus
	logl = mixLogl( all.data, all.parms, allPis, allMus, allLogDens, allLogls);

	//bundling up things to return
	//first the fitted pis
	double *tmpPi = REAL( Rpis);
	for( int i=0; i<all.data.nObs; i++)
		for( int g=0; g<all.data.nG; g++)
			tmpPi[MATREF(i,g,all.data.nObs)] = allPis.at(i).at(g);
	//the fitted expectations
	double *tmpMus = REAL( Rmus);
	for( size_t i=0; i<allMus.size(); i++)
		tmpMus[i] = allMus.at(i);
	//the log conditional densities
	double *tmpDens = REAL( RlogDens);
	for( int g=0; g<all.data.nG;g++)
		for( int i=0; i<all.data.nObs; i++)
				tmpDens[MATREF(i,g,all.data.nObs)] = allLogDens.at(i).at(g);
	//the logl contributions
	double *tmplogls = REAL( Rlogli);
	for( int i=0; i<all.data.nObs; i++)
		tmplogls[i] = allLogls.at(i);
	//the logl
	SEXP Rres;	//R object to return -- it is the logl!
	Rres = allocVector(REALSXP,1);
	double *R_pointer = REAL(Rres);
	R_pointer[0] = logl;

	return( Rres);
}
}

double invLogit( double x)
{
	double tmp;
	tmp = exp( x);
	tmp = tmp / (1+tmp);
	return( tmp);
}

double mixLogl( const myData &dat, const myParms &parms, vector< vector<double> > &allPis, vector<double> &allMus, vector< vector<double> > &allLogDens, vector<double> &allLogls)
{
	vector<double> logPis(dat.nG, dat.NAnum);//, pis( dat.nG, dat.NAnum);
	double res=0.0, pen=0.0;
	double wi, logli, peni;	//for coding purposes really
	vector<double> wij( dat.nG, dat.NAnum);	//for coding purposes
	int m;	//for coding purposes

	//calculate fitted values (constant over i)
	calcMuFits( allMus, dat, parms);
	for( int i=0; i<dat.nObs; i++){
		calcLogPis( logPis, allPis.at(i), dat, parms, i);
		calcLogCondDens( allLogDens.at(i), allMus, dat, parms, i);
		logli = calcMixSum( logPis, allLogDens.at(i), wi, wij, m);
		res += logli;
		peni = calcPiPen( logPis, dat, parms);
		pen += peni;
		allLogls.at(i) = logli + peni;
	}
//	penTau = calcTauPen( dat, parms);

	res += pen;
//	res += penTau;

	return( res);
}

double calcTauPen( const myData &dat, const myParms &parms)
{
	double penTau = 0.0;
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);

	parms.getAllTaus( newTau, dat);

	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++)
			penTau += -newTau.at(MATREF(g,s,dat.nG))*newTau.at(MATREF(g,s,dat.nG)) / (2*parms.sd*parms.sd);
	return( penTau);
}

double calcPiPen( const vector<double> &logPis, const myData &dat, const myParms &parms)
{
	double pen=0.0;

	for( int g=0; g<dat.nG; g++)
		pen += logPis.at(g);
	pen *= parms.conc;

	return( pen);
}

void calcLogPis( vector<double> &logPis, vector<double> &pis, const myData &dat, const myParms &parms, int i)
{
	vector<double> lp((dat.nG-1),0.0);
	double sumlp=0.0, sumpi=0.0;

	lp.assign( (dat.nG-1), 0.0);
	for( int k=0; k<(dat.nG-1); k++){
		for( int p=0; p<dat.np; p++)
			lp.at(k) += parms.Beta[MATREF(k,p,(dat.nG-1))] * dat.X[MATREF(i,p,dat.nObs)];
		lp.at(k) = exp( lp.at(k));
		sumlp += lp.at(k);
	}
	for( int k=0; k<(dat.nG-1); k++){
		pis.at(k) = lp.at(k) / ( 1+sumlp);
		sumpi += pis.at(k);
	}
	pis.at(dat.nG-1) = 1-sumpi;
	for( int k=0; k<dat.nG; k++)
		logPis.at(k) = log( pis.at(k));

}

void calcLogCondDens( vector<double> &condDens, const vector<double> &fits, const myData &dat, const myParms &parms, int i)
{
	vector<double> condDensSG( dat.nG*dat.nS, dat.NAnum);

	//calcualte the G*S log conditional densities
	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++)
			condDensSG.at(MATREF(g,s,dat.nG)) = logBernoulli( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF(g,s,dat.nG)));
	//calculate the G log conditional densities (under independence)
	for( int g=0; g<dat.nG; g++){
		condDens.at(g) = 0.0;
		for( int s=0; s<dat.nS; s++)
			condDens.at(g) += condDensSG.at(MATREF(g,s,dat.nG));
	}
}

void calcMuFits( vector<double> &fits, const myData &dat, const myParms &parms)
{
	//fits is a G*S matrix of the fitted values
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);
	vector<double> lps(dat.nG*dat.nS, dat.NAnum);

	//calculate sum-to-zero taus
	parms.getAllTaus( newTau, dat);
	//calcualte the G*S fits
	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++){
			lps.at(MATREF(g,s,dat.nG)) = parms.Alpha[s] + newTau[MATREF(g,s,dat.nG)];
			fits.at(MATREF(g,s,dat.nG)) = invLogit( lps.at(MATREF(g,s,dat.nG)));
		}
}

double logBernoulli( double y, double mu)
{
	double tmp;
	if( y==1){
		tmp = log( mu);
		return( tmp);
	}
	tmp = log( 1-mu);
	return( tmp);
}

double calcMixSum( const vector<double> &logPis, const vector<double> &condDens, double &wi, vector<double> &wij, int &maxg)
{
//function used in both mixLogl and loglDerivs
	vector<double> summands(logPis.size(), 0.0);
	double max, res=0.0;

	max = logPis.at(0) + condDens.at(0);
	maxg = 0;
	for( size_t i=0; i<logPis.size(); i++){
		summands.at(i) = logPis.at(i) + condDens.at(i);
		if( summands.at(i) > max){
			max = summands.at(i);
			maxg = i;
		}
	}
	wi = 0.0;
	for( size_t g=0; g<summands.size(); g++){
		wij.at(g) = exp( summands.at(g) - max);
		wi += wij.at(g);
	}
	res = log( wi);
	res += max;

	return( res);
}

void loglDerivs( const myData &dat, const myParms &parms, myDerivs &derivs)
{
	vector<double> logPis(dat.nG, dat.NAnum), pis( dat.nG, dat.NAnum);
	vector<double> logCondDens( dat.nG, dat.NAnum);
	double wi;
	vector<double> wij( dat.nG, dat.NAnum);
	vector<double> fits(dat.nG*dat.nS, dat.NAnum);
	int m;	//location of the maximum group contribution
	double tmp;

	vector<double> muDerivsI( dat.nG*dat.nS, dat.NAnum);
	vector<double> etaDerivsI( dat.nG*dat.nS, dat.NAnum);
	vector<double> alphaDerivsI( dat.nS, dat.NAnum);
	vector<double> tauDerivsI( (dat.nG-1)*dat.nS, dat.NAnum);
	vector<double> piDerivsI( dat.nG, dat.NAnum);
	vector<double> betaDerivsI( (dat.nG-1)*dat.np, dat.NAnum);

	vector<double> tmpPiDerivs( dat.nG*dat.nS, 0.0);

	calcMuFits( fits, dat, parms);
	derivs.zeroDerivs( dat);
	for( int i=0; i<dat.nObs; i++){
		//calculating the w_i and the {w_ig}
		calcLogPis( logPis, pis, dat, parms, i);
		calcLogCondDens( logCondDens, fits, dat, parms, i);
		tmp = calcMixSum( logPis, logCondDens, wi, wij, m);
		//calc deriv w.r.t. mu and the eta (all of them)
		calcDerivMu( muDerivsI, fits, dat, parms, wi, wij, m, i);
		calcDerivEtaMu( etaDerivsI, dat, muDerivsI, fits);
		//calc deriv w.r.t. alpha and then tau
		calcAlphaDeriv( alphaDerivsI, etaDerivsI, dat);
		calcTauDeriv( tauDerivsI, etaDerivsI, dat, parms);
		//calc deriv w.r.t. beta
		calcPiDeriv( piDerivsI, dat, parms, pis, wi, wij, m);
		calcBetaDeriv( betaDerivsI, piDerivsI, pis, dat, i);

		derivs.updateDerivs( dat, alphaDerivsI, tauDerivsI, betaDerivsI);
		derivs.updateScores( dat, alphaDerivsI, tauDerivsI, betaDerivsI, i);
	}
}

void calcDerivMu( vector<double> &muDerivs, const vector<double> &fits, const myData &dat, const myParms &parms, const double wi, const vector<double> &wij, const int m, const int i)
{
	//muDerivs is a GxS matrix of first derivatives
	vector<double> tmpDerivs( dat.nG*dat.nS, 0.0);

	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++)
			tmpDerivs.at(MATREF(g,s,dat.nG)) = logBernDer( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF(g,s,dat.nG)));

	for( int s=0; s<dat.nS; s++){
		muDerivs.at(MATREF(m,s,dat.nG)) = 0.0;
		for( int g=0; g<dat.nG; g++){
			if( m!=g){
				muDerivs.at(MATREF(g,s,dat.nG)) = tmpDerivs.at(MATREF(g,s,dat.nG)) * wij.at(g) / wi;
				muDerivs.at(MATREF(m,s,dat.nG)) -= tmpDerivs.at(MATREF(m,s,dat.nG)) * wij.at(g) / wi;
			}
			else
				muDerivs.at(MATREF(g,s,dat.nG)) += tmpDerivs.at(MATREF(g,s,dat.nG));
		}
	}
}

void calcDerivEtaMu( vector<double> &etaDerivsI, const myData &dat, const vector<double> &muDerivsI, const vector<double> &fits)
{
	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++)
			etaDerivsI.at(MATREF(g,s,dat.nG)) = fits.at(MATREF(g,s,dat.nG)) * (1-fits.at(MATREF(g,s,dat.nG))) * muDerivsI.at(MATREF(g,s,dat.nG));
}

double logBernDer( double y, double mu)
{
	double tmp;
	if( y==1){
		tmp = 1/mu;
		return( tmp);
	}
	if( y==0){
		tmp = -1/(1-mu);
		return( tmp);
	}
	return( log( -1));	//to give an error
}

void calcAlphaDeriv( vector<double> &alphaDerivsI, const vector<double> &etaDerivs, const myData &dat)
{
	alphaDerivsI.assign(alphaDerivsI.size(), 0.0);
	for( int s=0; s<dat.nS; s++)
		for( int g=0; g<dat.nG; g++)
			alphaDerivsI.at(s) += etaDerivs.at(MATREF(g,s,dat.nG));
}

void calcTauDeriv( vector<double> &tauDerivsI, const vector<double> &etaDerivs, const myData &dat, const myParms &parms)
{
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);

	tauDerivsI.assign(tauDerivsI.size(), 0.0);
	for( int s=0; s<dat.nS; s++){
		for( int g=0; g<(dat.nG-1); g++){
			tauDerivsI.at(MATREF(g,s,(dat.nG-1))) = etaDerivs.at(MATREF(g,s,dat.nG));
			tauDerivsI.at(MATREF(g,s,(dat.nG-1))) -= etaDerivs.at(MATREF((dat.nG-1),s,dat.nG));
		}
	}
}

void calcTauPenDeriv( vector<double> &tauDerivsI, const myData &dat, const myParms &parms)
{
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);

	tauDerivsI.assign(tauDerivsI.size(), 0.0);
	parms.getAllTaus( newTau, dat);
	for( int s=0; s<dat.nS; s++)
		for( int g=0; g<(dat.nG-1); g++){
			tauDerivsI.at(MATREF(g,s,(dat.nG-1))) += -newTau.at( MATREF(g,s,dat.nG)) / (parms.sd*parms.sd);
			tauDerivsI.at(MATREF(g,s,(dat.nG-1))) += newTau.at( MATREF((dat.nG-1),s,dat.nG)) / (parms.sd*parms.sd);
	}
}

void calcPiDeriv( vector<double> &piDerivsI, const myData &dat, const myParms &parms, const vector<double> &pis, const double wi, const vector<double> &wig, int m)
{
	vector<double> wigDerivs(dat.nG, 0.0);

	for( int g=0; g<dat.nG; g++){
		if( g!=m)
			piDerivsI.at(g) = wig.at(g) / (wi*pis.at(g));
	}
	piDerivsI.at(m) = 1 / pis.at(m);
	for( int g=0; g<dat.nG; g++)
		if( g!=m)
			piDerivsI.at(m) -= wig.at(g) / (wi*pis.at(m));

	for( int g=0; g<dat.nG; g++)
		piDerivsI.at(g) += parms.conc / pis.at(g);
}

void calcBetaDeriv( vector<double> &betaDerivsI, const vector<double> &piDerivsI, const vector<double> &pis, const myData &dat, int i)
{
	vector<double> dpideta( (dat.nG*(dat.nG-1)), 0.0);
	vector<double> dldeta( (dat.nG-1), 0.0);

	//derivs of pi w.r.t. eta (all (G-1) of 'em)
	for( int g=0; g<(dat.nG-1); g++){
		dpideta.at(MATREF(g,g,dat.nG)) += pis.at(g);
		for( int h=0; h<(dat.nG-1); h++)
			dpideta.at(MATREF(g,h,dat.nG)) += -pis.at(g) * pis.at(h);
	}
	for( int g=0; g<(dat.nG-1); g++){
		dpideta.at(MATREF((dat.nG-1),g,dat.nG)) = 0.0;
		for( int h=0; h<(dat.nG-1); h++)
			dpideta.at(MATREF((dat.nG-1),g,dat.nG)) -= dpideta.at(MATREF(h,g,dat.nG));
	}

	//deriv is dldpi X dpideta X detadbeta_h, for lp number h (of course)
	//logl_i w.r.t eta first: a 1x(G-1) vector
	for( int h=0; h<(dat.nG-1); h++)
		for( int g=0; g<dat.nG; g++)
			dldeta.at(h) += piDerivsI.at(g)*dpideta.at(MATREF(g,h,dat.nG));

	//now for each of the beta_h vectors
	betaDerivsI.assign( betaDerivsI.size(), 0.0);
	for( int h=0; h<(dat.nG-1); h++)
		for( int p=0; p<dat.np; p++)
			betaDerivsI.at(MATREF(h,p,(dat.nG-1))) += dldeta.at(h)*dat.X[MATREF(i,p,dat.nObs)];
}

double GSoptimise( allClasses &all)
{
	double *vmminGrad, *vmminParms, *oldParms;	//arrays to pass to vmmin
	vmminParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	oldParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	vmminGrad = (double *) R_alloc(all.parms.nTot,sizeof(double));
	int *myMask;
	vector<int> vecMask(all.parms.nTot, 0);
	double vmminLogl[1];

//	cout << "Gauss-Seidel iterations (" << all.contr.maxitGS << " of them)" << endl;
	Rprintf( "Gauss-Seidel iterations (%u of them) \n", all.contr.maxitGS);
	int kount=0, offs=0;
	do{
		all.parms.getArray( oldParms, all.data);
		Rprintf( "Habitat parameters...");
		for( int g=0; g<(all.data.nG-1); g++){
			//set mask
			offs = all.data.nS + (all.data.nG-1)*all.data.nS;
			for( int p=0; p<all.data.np; p++)
				vecMask.at(offs+MATREF(g,p,(all.data.nG-1))) = 1;
			myMask = &vecMask[0];
			//optimise for beta_g
			all.parms.getArray( vmminParms, all.data);
			vmmin( all.parms.nTot, vmminParms, vmminLogl, optimise_function, gradient_function, all.contr.maxitInner, all.contr.traceGS, myMask, all.contr.abstol, all.contr.reltolInner, all.contr.nReport, &all, &all.contr.fnKount, &all.contr.grKount, &all.contr.ifail);
			//update parameters
			all.parms.update( vmminParms, all.data);
			//tidy up
			vecMask.assign( vecMask.size(),0);
		}
		Rprintf( "Species parameters...");
		for( int s=0; s<all.data.nS; s++){
			//set mask
			vecMask.at(s) = 1;
			for( int g=0; g<(all.data.nG-1); g++)
				vecMask.at(all.data.nS+MATREF(g,s,(all.data.nG-1))) = 1;
			myMask = &vecMask[0];
			//optimise for alpha_s and tau_s
			all.parms.getArray( vmminParms, all.data);
			vmmin( all.parms.nTot, vmminParms, vmminLogl, optimise_function, gradient_function, all.contr.maxitInner, all.contr.traceGS, myMask, all.contr.abstol, all.contr.reltolInner, all.contr.nReport, &all, &all.contr.fnKount, &all.contr.grKount, &all.contr.ifail);
			//update parameters
			all.parms.update( vmminParms, all.data);
			//tidy up
			vecMask.assign( vecMask.size(),0);
		}
		kount++;
		all.parms.getArray( vmminParms, all.data);
//		cout.precision(5);
//		cout << "done!\n" << "GS iter: " << kount << " logl: " << fixed << vmminLogl[0] << "\n";
		Rprintf( "done! \n GS iter: %u \t logl: %10.5f \n", kount, vmminLogl[0]);
	}while( (!converged(oldParms, vmminParms, all.contr, all.parms.nTot)) && (kount < all.contr.maxitGS));//all.contr.maxitOuter));
	all.parms.update( vmminParms, all.data);
	gradient_function(all.parms.nTot,vmminParms, vmminGrad, &all);
	all.derivs.update( vmminGrad, all.data);
	if( kount < all.contr.maxitGS)
		*all.contr.conv = 1;
	else
		*all.contr.conv = 0;

	return( vmminLogl[0]);
}

double ALLoptimise( allClasses &all)
{
	double *vmminGrad, *vmminParms, *oldParms;	//arrays to pass to vmmin
	vmminParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	oldParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	vmminGrad = (double *) R_alloc(all.parms.nTot,sizeof(double));
	int *myMask;
	vector<int> vecMask(all.parms.nTot, 1);
	double vmminLogl[1];

//	cout << "Quasi-Newton iterations\n";
	Rprintf( "Quasi-Newton iterations\n");
	all.parms.getArray( oldParms, all.data);
	myMask = &vecMask[0];
	//optimise
	all.parms.getArray( vmminParms, all.data);
	vmmin( all.parms.nTot, vmminParms, vmminLogl, optimise_function, gradient_function, all.contr.maxitQN, all.contr.traceQN, myMask, all.contr.abstol, all.contr.reltolInner, all.contr.nReport, &all, &all.contr.fnKount, &all.contr.grKount, &all.contr.ifail);
	//update parameters
	all.parms.update( vmminParms, all.data);
	gradient_function(all.parms.nTot,vmminParms, vmminGrad, &all);
	all.derivs.update( vmminGrad, all.data);

	return( vmminLogl[0]);
}

bool converged( double *oldP, double *newP, const myOptContr &contr, int nTot)
{
	double tmp, eps=1e-5;

	for( int i=0; i<nTot; i++){
		tmp = fabs(newP[i]-oldP[i]);
		tmp /= fabs(oldP[i])+eps;
		if( tmp > contr.reltol)
			return( false);
	}
	return( true);

}

double optimise_function(int n, double *par, void *ex)
{
	allClasses *all = (allClasses *) ex;
	double logl;
	vector< vector<double> > allPis( all->data.nObs, vector<double> (all->data.nG, all->data.NAnum));	//(nObs x nG) matrix --2D
	vector<double> allMus( all->data.nG*all->data.nS, all->data.NAnum);
	vector< vector<double> > allLogDens( all->data.nObs, vector<double> (all->data.nG, all->data.NAnum));
	vector<double> allLogl( all->data.nObs, all->data.NAnum);

	all->parms.update( par, all->data);
	logl = mixLogl( all->data, all->parms, allPis, allMus, allLogDens, allLogl);
	return( (0.0-logl));
}

void gradient_function(int n, double *par, double *gr, void *ex)
{
	allClasses *all = (allClasses *) ex;

	loglDerivs( all->data, all->parms, all->derivs);

	all->derivs.getArray(gr, all->data);

	for( int i=0; i<n; i++)
		gr[i] = 0-gr[i];
}
