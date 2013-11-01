\name{simRCPdata}
\alias{simRCPdata}
\title{Simulates from a regimix model}
\description{Simulates a data set from a mixture-of-experts model for RCP (for region of common profile) types.
}
\usage{ simRCPdata(H=3, S=20, p=3, n=200, alpha=NULL, tau=NULL, beta=NULL, X=NULL)
}
\arguments{
\item{ H}{Integer giving the number of RCPs}
\item{ S}{Integer giving the number of species}
\item{ p}{Integer giving the number of covariates (including the intercept)}
\item{ n}{Integer giving the number of observations (sites)}
\item{ alpha}{Numeric vector of length S. Specifies the mean prevalence for each species, on the logit scale.}
\item{ tau}{Numeric matrix of dimension c(H-1,S). Specifies each species difference from the mean to each RCPs mean for the first H-1 RCPs. The last RCP means are calculated using the sum-to-zero constraints}
\item{ beta}{Numeric matrix of dimension c(H-1,p). Specifies the RCP's dependence on the covariates (in X)}
\item{ X}{Numeric matrix of dimension c(n,p). Specifies the covariates for the RCP model. Must include the intercept, if one is wanted}
}
\value{
A data frame that contains the outcomes (species data) and the covariates (environmental data). This data.frame has a number of special attirbutes, which are information about the model underlying the data. They are:
\item{ RCPs}{the true, but unobserved, RCP types}
\item{ pis}{the true prior probabilities}
\item{ alpha}{the species overall prevalences, on linear predictor scale}
\item{ tau}{the deviation from alpha for each RCP type, on linear predictor scale}
\item{ beta}{the parameters controlling how the RCP types depend on the covariates}
\item{ mu}{the probabilities of each species occuring in each RCP type}
}
\author{Scott D. Foster}
\references{Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics.
}
\examples{
#generates synthetic data in Foster et al (2012)
set.seed( 151)
n <- 150
S <- 100
H <- 3
X <- as.data.frame( cbind( x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10)))
pols <- list()
pols[[1]] <- poly( X$x1, degree=3)
pols[[2]] <- poly( X$x2, degree=3)
X <- as.matrix( cbind( 1, X, pols[[1]], pols[[2]]))
colnames( X) <- c("const", 'x1', 'x2', paste( "x1",1:3,sep='.'), paste( "x2",1:3,sep='.'))
p <- ncol( X[,-(2:3)])
alpha <- rnorm( S)
tau.var <- 0.5
b <- sqrt( tau.var/2)
tau <- matrix( rexp( n=(H-1)*S, rate=1/b) - rexp( n=(H-1)*S, rate=1/b), nrow=H-1, ncol=S)
beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0, -14.2, -28.3,
  -36.8, -17.8, -92.9,-2.7), nrow=H-1, ncol=p)
set.seed(121)
simDat <- simRCPdata( H=H, S=S, p=p, n=n, alpha=alpha, tau=tau, beta=beta, X=X[,-(2:3)])
}
\keyword{misc}
