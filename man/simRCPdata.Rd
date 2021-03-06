\name{simRCPdata}
\alias{simRCPdata}
\title{Simulates from a regimix model}
\description{Simulates a data set from a mixture-of-experts model for RCP (for region of common profile) types.
}
\usage{
 simRCPdata(nRCP=3, S=20, n=200, p.x=3, p.w=0, alpha=NULL, tau=NULL, beta=NULL,
                        gamma=NULL, logDisps=NULL, powers=NULL, X=NULL, W=NULL,
                                                  offset=NULL, dist="Bernoulli")
}
\arguments{
\item{ nRCP}{Integer giving the number of RCPs}
\item{ S}{Integer giving the number of species}
\item{ n}{Integer giving the number of observations (sites)}
\item{ p.x}{Integer giving the number of covariates (including the intercept) for the model for the latent RCP types}
\item{ p.w}{Integer giving the number of covariates (excluding the intercept) for the model for the species data}
\item{ alpha}{Numeric vector of length S. Specifies the mean prevalence for each species, on the logit scale}
\item{ tau}{Numeric matrix of dimension c(nRCP-1,S). Specifies each species difference from the mean to each RCPs mean for the first nRCP-1 RCPs. The last RCP means are calculated using the sum-to-zero constraints}
\item{ beta}{Numeric matrix of dimension c(nRCP-1,p.x). Specifies the RCP's dependence on the covariates (in X)}
\item{ gamma}{Numeric matrix of dimension c(n,p.w). Specifies the species' dependence on the covariates (in W)}
\item{ logDisps}{Logartihm of the (over-)dispersion parameters for each species for negative binomial, Tweedie and Normal models}
\item{ powers}{Power parameters for each species for Tweedie model}
\item{ X}{Numeric matrix of dimension c(n,p.x). Specifies the covariates for the RCP model. Must include the intercept, if one is wanted. Default is random numbers in a matrix of the right size.}
\item{ W}{Numeric matrix of dimension c(n,p.w). Specifies the covariates for the species model. Must \emph{not} include the intercept. Unless you want it included twice. Default is to give random levels of a two-level factor.}
\item{ offset}{Numeric vector of size n. Specifies any offset to be included into the species level model.}
\item{ dist}{Text string. Specifies the distribution of the species data. Current options are "Bernoulli" (default), "Poisson", "NegBin", "Tweedie" and "Normal".}
}
\value{
A data frame that contains the outcomes (species data) and the covariates (environmental data and species-level covariates). This data.frame has a number of special attirbutes, which are information about the model underlying the data. They are:
\item{ RCPs}{the true, but unobserved, RCP types}
\item{ pis}{the true prior probabilities}
\item{ alpha}{the species overall prevalences, on linear predictor scale}
\item{ tau}{the deviation from alpha for each RCP type, on linear predictor scale}
\item{ beta}{the parameters controlling how the RCP types depend on the covariates}
\item{ gamma}{the parameters controlling how each species depends on the species-level covariates}
\item{ logDisps}{the logarithm of the dispersion parameter for each species}
\item{ mu}{the probabilities of each species occuring in each RCP type}
}
\author{Scott D. Foster}
\references{Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics.
}
\examples{
#generates synthetic data
set.seed( 151)
n <- 100
S <- 10
nRCP <- 3
my.dist <- "NegBin"
X <- as.data.frame( cbind( x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10)))
Offy <- log( runif( n, min=30, max=60))
pols <- list()
pols[[1]] <- poly( X$x1, degree=3)
#important to scale covariates so that regimix can get half-way decent starting values
pols[[2]] <- poly( X$x2, degree=3)
X <- as.matrix( cbind( 1, X, pols[[1]], pols[[2]]))
colnames( X) <- c("const", 'x1', 'x2', paste( "x1",1:3,sep='.'), paste( "x2",1:3,sep='.'))
p.x <- ncol( X[,-(2:3)])
p.w <- 3
W <- matrix(sample( c(0,1), size=(n*p.w), replace=TRUE), nrow=n, ncol=p.w)
colnames( W) <- paste( "w",1:3,sep=".")
alpha <- rnorm( S)
tau.var <- 0.5
b <- sqrt( tau.var/2)
#a double exponential for RCP effects
tau <- matrix( rexp( n=(nRCP-1)*S, rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b), nrow=nRCP-1, ncol=S)
beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0, -14.2, -28.3,
  -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)
logDisp <- log( rexp( S, 1))
set.seed(121)
simDat <- simRCPdata( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n, alpha=alpha, tau=tau,
  beta=beta, gamma=gamma, X=X[,-(2:3)], W=W, dist=my.dist, logDisp=logDisp, offset=Offy)
}
\keyword{misc}
