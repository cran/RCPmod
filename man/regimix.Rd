\name{regimix}
\alias{regimix}
\alias{regimix.multifit}
\title{Fits a regimix model.}
\description{ Fits a mixture-of-experts model to identify regions of similar community composition.
}
\usage{
 regimix( form.RCP = NULL, form.spp = NULL, data, nRCP = 3, dist="Bernoulli",
    offset=NULL, weights=NULL, control = list(), inits = "random2",
    titbits = TRUE, power=1.6)
 regimix.multifit( form.RCP = NULL, form.spp = NULL, data, nRCP = 3, dist="Bernoulli",
    offset=NULL, weights=NULL, control = list(), inits = "random2",
    titbits = FALSE, power=1.6, nstart=10, mc.cores=1)
}
\arguments{
\item{ form.RCP}{an object of class "formula" (or an object that can be coerced to that class). The left hand side of this formula specifies the columns of the data argument that are the species binary data. The right hind side of this formula specifies the dependence of the \code{region of common profile} (RCP) probabilities on covariates. An example formula is cbind(spp1,spp2,spp3)~1+poly(temp,3) where spp1, spp2 and spp3 are species labels (columns of data), and the RCP probabilities depend on a (link-)cubic polynomial of temp.}
\item{ form.spp}{an object of class "formula" (or an object that can be coerced to that class). The left hand side of this formula should be left empty (it is removed if it is not empty). The right hand side of this formula specifies the dependence of the species"'" data on covariates (typically different covariates to form.RCP to avoid confusing confounding). An example formula is ~gearType + timeOfDay, where gearType describes the different sampling gears and timeOfDay describes the time of the sample. Any intercept term in this model is removed (as this is dealt with in the species-specific RCP-specific intercepts).}
\item{ data}{an object of class "data.frame" (or one that can be coerced to that class). This data.frame has to contain all the data for the terms in the form.RCP, form.spp and offset arguments.}
\item{ nRCP}{an integer (or something that can be coerced to an integer). This argument specifies the number of RCPs that will be fitted.}
\item{ dist}{a character string describing the distribution of the species data. Must be one of "Bernoulli" (for binary data), "Poisson" and "NegBin" (for count data), "Tweedie" (for biomass, non-negative, data), and "Normal" (for quantity data). Note that this also specifies the link function used for the data. The links are:logit (for "Bernoulli"), log (for "Poisson", "NegBin" and "Tweedie") and identity (for "Normal"). Note that the computation strategy for the Tweedie model is slightly different than for the other models. The Tweedie models can also take a while to fit, even for small data sets, as calculating Tweedie densities can be quite labourious.}
\item{ offset}{a numeric vector of length nrow( data) that is included into the model as an offset. It is included into the conditional part of the model where conditioning is performed on the unobserved RCP type. Note that offsets cannot be included as part of the form.RCP or form.spp arguments -- only through this argument.}
\item{ weights}{a numeric vector of length nrow( data) that is used as weights in the log-likelihood calculations. If NULL (default) then all weights are assumed to be identically 1.}
\item{ control}{a list of control parameters for optimisation and calculation. See details.}
\item{ inits}{either a character string or a numeric vector. If a character string ("random2", "random", "hclust" or "noPreClust") then it gives the method to generate starting values. If a numeric vector then it specifies the values of alpha, tau, beta, gamma, and log( dispersions), in that order. It will be used unchecked if not a character. The default is "random2" which is described in Foster et al (in prep.). Other choices are "random", which is described in Foster et al (2013), "hclust" which is the same as "random2" but with no random component (also described in Foster et al (2013)). The "noPreClust" choice is designed for situations where running hclust() on the data is not feasible (due to data size) or not wanted -- instead the initial (design) matrix for assigning sites to RCPs is obtained by samples from a symmetric Dirichlet distribution with shape parameter 5. It is strongly advised that multiple starts with multiple (random) start locations are performed. The reason for this is that the log-likelihood surface can be fairly "bumpy" with multiple local maxima. Multiple starts guards somewhat against making inference from these local maximia. For regimix.multifit, inits should be the either "random" or "random2" (default).}
\item{ titbits}{either a boolean or a vector of characters. If TRUE (default for regimix(qv)), then some objects used in the estimation of the model"'"s parameters are returned in a list entitled "titbits" in the model object. Some functions, for example plot.regimix(qv) and predict.regimix(qv), will require some or all of these pieces of information. If titbits=FALSE (default for regimix.multifit(qv)), then an empty list is returned. If a character vector, then just those objects are returned. Possible values are:"Y" for the outcome matrix, "X" for the model matrix for the RCP model, "W" for the model matrix for the species-specific model, "offset" for the offset in the model, "wts" for the model weights, "form.RCP" for the formula for the RCPs, "form.spp" for the formula for the species-specific model, "control" for the control arguments used in model fitting, "dist" for the conditional distribution of the species data, and "power" for the power parameters used (only used in Tweedie models). Care needs to be taken when using titbits=TRUE in regimix.multifit(qv) calls as titbits is created for EACH OF THE MODEL FITS. If the data is large or if nstart is large, then setting titbits=TRUE may give users problems with memory. It is more efficient, from a memory perspective, to refit the "best" model using regimix(qv) after identifying it with regimix.multifit(qv). See examples for illustration about how to do this.}
\item{ power}{a numeric vector (length either 1 or the number of species) defining the power parameter to use in the Tweedie models. If length(power)==1, then the same power parameter is used for all species. If length(power)==No_species, then each species gets its own power parameter. Power values must be between 1 and 2, for computational reasons they should be well away from the boundary. The default is 1.6 as this has proved to be a good ball-park value for the fisheries data that the developer has previously analysed.}
\item{ nstart}{for regimix.multifit only. The number of random starts to perform for re-fitting. Default is 10, which will need increasing for serious use.}
\item{ mc.cores}{for regimix.multifit only. The number of cores to spread the re-fitting over.}
}
\details{
 A typical formula for use in the form.RCP argument will have the form (for example) cbind(spp1,spp2,spp3,spp4)~1+cov1+cov2*cov3. This signifies that there are 4 species to be used for RCP modelling and that the RCP types are dependent on cov1+cov2+cov3+cov2:cov3. See ?glm for a description of how the right hand side of the formula is expanded.

 Likewise a typical formula for use in the form.spp argument will have the form (for example) ~1+fac1+cov1. This signifies that the catchabilities of each species depends upon the levels of the factor fac1 and the covariate cov1. See ?glm for a description of how the right hand side of the formula is expanded.

 The computation strategy for the default method, which has been demonstrated to work for all data sets encountered thus far, is fully described in Foster et al (2013) and Foster et al (in prep.). We note however, that it is a good idea to standardise covariates prior to calling regimix. This is not formally required by the model, but the method to obtain starting values is based on a penalised likelihood, which is dependent on covariate scale. To recap, the starting values are dependent on covariate scale but the final solution (if it can be found given starting values) should not be.

 We do not, on purpose, provide residuals as a routine part of the model. Users should use the residuals.regimix(qv) function to obtain them. We do this as the type of residual needs to be specified (although we recommend type=="RQR" for routine use).

 Control arguments for optimisation generally follow those in optim(qv), although a few differences occur (e.g. "loglOnly"). The elements of the control list are

\describe{
\item{maxit}{The maximum number of iterations. Default is 500.}
\item{quiet}{Should any reporting be performed? Default is FALSE, for reporting. For regimix.multifit(), this indicates if the progress should be printed.}
\item{trace}{Non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information.}
\item{nreport}{The frequency of reports for optimisation. Default is 10 -- a report for 10th iteration.}
\item{reltol}{Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8.}
\item{optimise}{Should optimisation for estimation occur? If TRUE (default) optimisation will occur. If FALSE no optimisation is performed.}
\item{loglOnly}{Should the log-likelihood be caulcated? If TRUE (default) then log-likelihood is calculated and returned. If FALSE then the log-likelihood is not calculated for return.}
\item{derivOnly}{Should the scores be evaluated at the (final) parameter values. If TRUE (default) then they are calculated. If FALSE then they are not calculated.}
\item{penalty}{A numeric scalar. This is the concentration for the Dirichlet-inspired penalty for the prior probabilities. Values less than zero will be set to the default (0.1). Large values give more penalisation than small ones.}
\item{penalty.tau}{A numeric scalar. This is the penalty for the tau parameters in the species model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
\item{penalty.gamma}{A numeric scalar. This is the penalty for the gamma parameters in the species model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
\item{penalty.disp}{a two element vector. These are combined to form the penalty for the dispersion parameters (if any). The dispersions are assumed to come from a log-normal distribution with log-mean penalty.disp[1] and log-standard-deviation penalty.disp[2]. Defaults to c(10,sqrt(10)), which gives shrinkage towards 1 (the mode of the penalty). Note that for Normal models, where the dispersion alone defines the variance, a strong penalty may be required to keep parameters estimable.}
}

 For calls to regimix.multifit(), titbits is set to FALSE-- so no excess memory is used. If users want this information, and there is good reason to want it, then a call to regimix() with starting values given as the \code{best} fit's estimates should be used.
}
\value{
regimix returns an object of class \code{regimix} and regimix.multfit returns a list of objects of class \code{regimix}. The regimix class has several methods: coef, plot, predict, residuals, summary, and vcov. The regimix object consists of a list with the following elements:
\item{ AIC}{Akaike an information criterion for the maximised model.}
\item{ BIC}{Bayesian information criterion for the maximised model.}
\item{ call}{the call to the function.}
\item{ coefs}{a list of three elements, one each for the estimates for the species prevalence (alpha), the deviations from alpha for the first (nRCP-1) RCP (tau), and the (nRCP-1) sets of RCP regression coefficents (beta).}
\item{ conv}{the convergence code from the maximisation procedure. See ?optim for an explanation.}
\item{ dist}{the character string identifying the distribution used for the model.}
\item{ logCondDens}{an nObs by nRCP matrix specifying the probability of observing each sites"'" data, given each of the RCP types.}
\item{ logl}{the maximised log likelihood.}
\item{ mus}{an array of size nRCP x S x nRCP where each element of the first dimension is the fitted value for all the species in all the RCP types.}
\item{ n}{the number of samples.}
\item{ names}{the names of the species, and the names of the covariates for X and W.}
\item{ nRCP}{the number of RCPs.}
\item{ pis}{an n x nRCP matrix with each column giving the prior probabilities for the corresponding RCP type. Rows sum to one.}
\item{ postProbs}{an n x nRCP matrix with each column giving the posterior probabilities for the corresponding RCP type. Rows sum to one (as each site is assumed to be from one of the RCP types).}
\item{ p.w}{the number of covariates used in the species-specific model.}
\item{ p.x}{the number of covariates used in the RCP model}
\item{ S}{the number of species.}
\item{ scores}{a list of three elements. Structure corresponds to coefs.}
\item{ start.vals}{the values used to start the estimation procedure.}
\item{ titbits}{(if requested using the titbit argument, see above) other pieces of information, useful to developers, that users should not typically need to concern themselves with. However, this information is used by methods for regimix objects.}
}
\author{Scott D. Foster}
\references{
Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245

Foster, S.D., Lyons, M. and Hill, N. (in prep.) Ecological Groupings of Sample Sites in the presence of sampling artefacts.
}
\examples{
#simulate data
example( simRCPdata)  #generates Negative Binomial data
#fit the model
my.form.RCP <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
    ')',sep=''),
       '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
my.form.spp <- ~w.1+w.2+w.3
fm <- regimix( form.RCP = my.form.RCP, form.spp=my.form.spp, data = simDat,
                  dist="NegBin", nRCP = 3, inits = "random2", offset=offset)
\dontrun{
#fit the model using multiple starting values
fm <- regimix.multifit( form.RCP = my.form.RCP, form.spp=my.form.spp, data = simDat,
  dist="NegBin", nRCP = 3, inits = "random2", offset=offset, nstart=10, titbits=FALSE,
  mc.cores=1)
#sometimes the model 'mis-fits' and one or more of the RCP groups has no sites associated
#with it.  These need to be removed (based on the colSums of the posterior probabilities)
postProbSums <- t( sapply( fm, function(x) colSums( x$postProbs)))
#Identify those models with zero posterior prob classes
allGoodUns <- apply( postProbSums, 1, function(x) all(x!=0))
#subset the fits
fm.clean <- fm[allGoodUns]
#choose the model with the lowest BIC
goodUn <- which.min( sapply( fm.clean, BIC))
#Using the 'best' model, use regimix(qv) again to additional model output needed for other
#functions (e.g. plot.regimix(qv), predict.regimix(qv) and regiboot(qv)). Note that the
#model is not estimated again (see control argument of the following regimix(qv) call.
fm.final <- regimix( form.RCP = my.form.RCP, form.spp=my.form.spp, data = simDat,
  dist="NegBin", nRCP = 3, inits = unlist( fm.clean[[goodUn]]$coef),
  control=list(optimise=FALSE), offset=offset)
}
}
\keyword{misc}
