\name{AIC.regimix}
\alias{AIC.regimix}
\title{Information criterion for a regimix model.}
\description{Returns information criterion for regimix models.
}
\section{Method}{ AIC( object, \dots{}, k=2)
}
\arguments{
\item{ object}{an object obtained from fitting a region of common profile mixture model. Such as that generated from a call to regimix(qv).}
\item{ k}{the coefficient for the penalty in the information criterion. k=2 signifies Akaikes information criterion, k=log(object$n) corresponds to the Bayesian information criterion. If NULL (default) the AIC is used.}
\item{ ...}{ignored}
}
\value{A numeric scalar giving the information criterion.
}
\author{Scott D. Foster}
\keyword{misc}
