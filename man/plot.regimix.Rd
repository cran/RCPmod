\name{plot.regimix}
\alias{plot.regimix}
\title{Plots a regimix object.}
\description{Plots the residuals, along with simulation envelopes, from a regimix object.
}
\section{Method}{ plot( x, ..., nsim=100, alpha.conf=c(0.90,0.95,0.99))
}
\arguments{
\item{ x}{an object obtained from fitting a RCP mixture model. Such as that generated from a call to regimix(qv).}
\item{ nsim}{the number of simulations to use for the confidence interval. Default is 100, serious usage is likely to require more.}
\item{ alpha.conf}{the confidence level(s) to use in the residual plots. Default is c(0.90,0.95,0.99)}
\item{ ...}{ignored}
}
\details{Produces a plot with two panels. The first is for point-wise confidence intervals and the second is for approximate global intervals. See Foster et al (2012) for details.
}
\author{Scott D. Foster}
\references{Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics.
}
\keyword{misc}
