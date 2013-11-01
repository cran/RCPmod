\name{extractAIC.regimix}
\alias{extractAIC.regimix}
\title{Extracts the AIC for a fitted refimix model.}
\description{Computes the generalised AIC for a regimix model.
}
\section{Method}{ extractAIC( fit, scale, k=2, \dots{})
}
\arguments{
\item{ fit}{an object obtained from fitting a region of common profile mixture model. Such as that generated from a call to regimix(qv).}
\item{ scale}{ignored}
\item{ k}{the coefficient for the penalty in the information criterion. k=2 corresponds to Akaikes information criterion, k=log(object$n) corresponds to the Bayesian information criterion. Default is k=2 (AIC).}
\item{ ...}{ignored}
}
\value{A two element numeric vector. First element is the number of parameters in the model. The second is the information criterion.
}
\author{Scott D. Foster}
\keyword{misc}
