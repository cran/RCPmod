\name{coef.regimix}
\alias{coef.regimix}
\title{A regimix objects coefficients.}
\description{Returns coefficicients from a regimix object.
}
\section{Method}{ coef( object, ...)
}
\arguments{
\item{ object}{an object obtained from fitting a regions of common profile mixture model. Such as that generated from a call to regimix(qv).}
\item{ ...}{ignored}
}
\value{Returns a list of four elements, one each for the estimates for the species prevalence (alpha), the deviations from alpha for the first (nRCP-1) regional profiles (tau), the (nRCP-1) sets of region regression coefficents (beta), the coefficients for the species specific model (if specified in the model call), and the log of the dispersion parameter estimates (for negative binomial, Tweedie and Normal models).
}
\author{Scott D. Foster}
\keyword{misc}
