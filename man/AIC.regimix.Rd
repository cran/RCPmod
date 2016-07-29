\name{AIC.regimix}
\alias{AIC.regimix}
\alias{BIC.regimix}
\title{Information criterion for a regimix model.}
\description{Returns information criterion for regimix models.
}
\section{Method}{
 AIC( object, ..., k=2)

 BIC( object, ...)
}
\arguments{
\item{ object}{an object obtained from fitting a region of common profile mixture model. Such as that generated from a call to regimix(qv).}
\item{ ...}{ignored}
\item{  k}{the coefficient for the penalty in the information criterion. k=2 signifies Akaikes information criterion, k=log(object$n) corresponds to the Bayesian information criterion. If NULL (default) the AIC is used.}
}
\value{A numeric scalar giving the Bayesian information criterion.
}
\author{Scott D. Foster}
\keyword{misc}
