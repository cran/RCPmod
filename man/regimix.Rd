\name{regimix}
\alias{regimix}
\title{Fits a regimix model.}
\description{ Fits a mixture-of-experts model to identify regions of similar community composition.
}
\usage{ regimix( form=NULL, data, nRCP=3, control=list(), inits="random")
}
\arguments{
\item{ form}{an object of class "formula" (or one that can be coerced to that class). The left hand side of this formula specifies the columns of the data argument that are the species binary data. The right hind side of this formula specifies the dependence of the \code{region of common profile} (RCP) probabilities on covariates. An example formula is cbind(spp1,spp2,spp3)~1+poly(temp,3) where spp1, spp2 and spp3 are species labels (columns of data), and the RCP probabilities depend on a (transformed) cubic polynomial of temp. No offset term can be accounted for, please do not enter one -- the behaviour is unknown and could cause your computer to have awkward moments.}
\item{ data}{an object of class "data.frame" (or one that can be coerced to that class). This data.frame contains all the data for the terms in the form argument}
\item{ nRCP}{an integer (or something that can be coerced to an integer). This argument specifies the number of RCPs that will be fitted.}
\item{ control}{a list of control parameters for optimisation and calculation. See details.}
\item{ inits}{either a character string or a numeric vector. If a character string ("hclust" or "random") then it gives the method to generate starting values. If a numeric vector then it specifies the values of alpha, tau, and beta (in that order). It will be used unchecked if not a character.}
}
\details{
 A typical formula for use in the form argument will have the form (for example) cbind(spp1,spp2,spp3,spp4)~1+cov1+cov2*cov3. This signifies that there are 4 species to be used for RCP modelling and that the RCP types are dependent on cov1+cov2+cov3+cov2:cov3. See ?glm for a description of how the right hand side of the formula is expanded.

 The computation strategy for the default method, which has been demonstrated to work, is fully described in Foster et al (2012).

 Control arguments for optimisation generally follow those in optim(qv), although a few differences occur (e.g. maxitInner). The elements of the control list are

\describe{
\item{maxit}{The maximum number of iterations. Default is 500.}
\item{maxGSiters}{The maximum number of Gauss-Seidel iterations in a hybrid Gauss-Seidel -- quasi-Newton optimisation scheme. This specifies the number of Gauss-Seidel iterations to use before switching to quasi-Newton. Default is 0.}
\item{maxitInner}{The maximum number of iterations for each update within a Gauss-Seidel update. Default is 100. This argument is redundant if maxGSiters=0.}
\item{trace}{Non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information.}
\item{nreport}{The frequency of reports for optimisation. Default is 1 -- a report for every iteration.}
\item{reltol}{Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8.}
\item{loglOnly}{Should only the log-likelihood be caulcated? If TRUE then optimisation does not occur. If FALSE (default) then optimisation is performed.}
\item{penalty}{A numeric scalar. This is the concentration for the Dirichlet-inspired penalty for the prior probabilities. Values less than zero will be set to the default, 0.01.}
}
}
\value{
regimix returns an object of class \code{regimix}. The regimix class has several methods: coef, plot, predict, residuals, summary, and vcov. The regimix object consists of a list with the following elements:
\item{ n}{the number of samples.}
\item{ S}{the number of species.}
\item{ p}{the number of covariates used to delineate the regitats.}
\item{ nRCP}{the number of RCPs}
\item{ logl}{the maximised log likelihood.}
\item{ coefs}{a list of three elements, one each for the estimates for the species prevalence (alpha), the deviations from alpha for the first (nRCP-1) RCP (tau), and the (nRCP-1) sets of RCP regression coefficents (beta).}
\item{ scores}{a list of three elements. Structure corresponds to coefs.}
\item{ pis}{an n x nRCP matrix with each column giving the prior probabilities for the corresponding RCP type. Rows sum to one.}
\item{ postProbs}{an n x nRCP matrix with each column giving the posterior probabilities for the corresponding RCP type. Rows sum to one.}
\item{ mus}{a matrix of size nRCP x S where each row is the probability of each species being observed in that RCP type.}
\item{ residuals}{deviance residuals for each site.}
\item{ BIC}{Bayesian information criterion for the maximised model.}
\item{ AIC}{Akaike an information criterion for the maximised model.}
\item{ start.vals}{the values used to start the estimation procedure.}
\item{ call}{the call to the function.}
\item{ titbits}{other pieces of information, useful to developers, that users should not need to concern themselves with.}
}
\author{Scott D. Foster}
\references{Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics.
}
\examples{
#simulate data
example( simRCPdata)
#fit the model
my.form <- paste( paste( paste(
  'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
    ')',sep=''),
       '~x1.1+x1.2+x1.3+x2.1+x2.2+x2.3',sep='')
my.cont <- list( maxGSiters = 0, maxit=1000, nreport=10, penalty=0.01)
fm <- regimix(form = my.form, data = simDat, nRCP = 3, control = my.cont, inits = "hclust")
}
\keyword{misc}
