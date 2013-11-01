\name{coef.regimix}
\alias{coef.regimix}
\title{A regimix objects coefficients.}
\description{Returns coefficicients from a regimix object.
}
\section{Method}{ coef( object, \dots{})
}
\arguments{
\item{ object}{an object obtained from fitting a regions of common profile mixture model. Such as that generated from a call to regimix(qv).}
\item{ ...}{ignored}
}
\value{Returns a list of three elements, one each for the estimates for the species prevalence (alpha), the deviations from alpha for the first (nRCP-1) regional profiles (tau), and the (nRCP-1) sets of region regression coefficents (beta).
}
\author{Scott D. Foster}
\keyword{misc}
