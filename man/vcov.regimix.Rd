\name{vcov.regimix}
\alias{vcov.regimix}
\title{Variance matrix for a regimix object.}
\description{Calculates variance-covariance matrix from a regimix object.
}
\section{Method}{ vcov( object, ..., method='simple', nboot=100)
}
\arguments{
\item{ object}{an object obtained from fitting a RCP (for region of common profile) mixture model. Such as that generated from a call to regimix(qv).}
\item{ method}{the method to calculate the variance-covariance matrix. Options are:'simple' (default) and \code{Richardson} from the numDeriv package (finite difference approximations, see jacobian(qv) from numDeriv package); \code{EmpInfo}; and \code{boot} for empirical bootstrapping. The \code{EmpInfo} and \code{boot} methods are purely experimental. In fact, they are unlikely to work. The \code{simple} and \code{Richardson} methods will work for the regression parameters and \bold{may} work for the species specific parameters.}
\item{ nboot}{the number of bootstrap samples to take for the bootstrap estimation. Argument is ignored if method != \code{boot}}
\item{ ...}{ignored}
}
\details{
 If method is either \code{simple} or \code{Richardson}, then the estimates variance matrix is based on a finite difference approximation to the observed information matrix. This method appears to be the most robust for the parameters relating the RCP types to the covariates.

 If method is \code{EmpInfo} then an empirical approximation to the Fisher Information is used as the basis for the variance matrix estimate. These variances are only available for the parameters that relate the RCP types to the covariates. This method is \bold{experimental}! Initial results are not good.

 If method is \code{boot} then an empirical bootstrap is used to estimate the variance matrix. This method is \bold{experimental} and initial results suggest that the resulting variances are downwardly biased.
}
\value{A square matrix of size equal to the number of parameters. It contains the variance matrix of the parameter estimates. If method=='EmpInfo', then only the sub-matrix corresponding to those parameters relating the RCP types to the covariates have non-NA elements.
}
\author{Scott D. Foster}
\keyword{misc}
