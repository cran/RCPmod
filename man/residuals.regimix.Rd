\name{residuals.regimix}
\alias{residuals.regimix}
\title{Residuals for a regimix object.}
\description{ Returns the (absolute) deviance residuals or the randomised quantile residuals from a regimix object.
}
\section{Method}{ residuals( object, ..., type="RQR", quiet=FALSE)
}
\arguments{
\item{ object}{an object obtained from fitting a RCP mixture model. Such as that generated from a call to regimix(qv).}
\item{ ...}{ignored}
\item{ type}{either "RQR" (the default) or "deviance" (the old default). See details.}
\item{ quiet}{should information be printed during the function? quiet=FALSE gives the information, whereas quiet=TRUE does not.}
}
\details{
 The randomised quantile residuals ("RQR", from Dunn and Smyth, 1996) are defined by their marginal distribution function (marginality is over other species observations within that site; see Foster et al, in prep). The result is one residual per species per site and they all should be standard normal variates. Within a site they are likely to be correlated (as they share a common latent factor), but across sampling locations they will be independent.

 The deviance residuals (as used here), are actually just square root of minus two times the log-likelihood contribution for each sampling location. We do not subtract the log-likelihood of the saturated model as, at the time of writing, we are unsure what this log-likelihood should be (latent factors confuse things here). This implies that the residuals will not have mean zero and their variance might also be heteroskedastic. This was not realised when writing the original RCP paper (Foster et al, 2013), obviously. We still believe that these residuals have some utility, but we are unsure where that utility stops. For general useage, the "RQR" residuals should probably be preferred.
}
\value{
\item{ For type=="RQR", a number-of-sites by number-of-species matrix with the randomised quantile residuals, which should be distributed as a standard normal variate.}{}
\item{ For type=="deviance" a numeric vector of size object$n containing the deviance residuals.}{}
}
\author{Scott D. Foster}
\references{
Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.

Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245

Foster, S.D., Lyons, M. and Hill, N. (in prep.) Ecological Groupings of Sample Sites in the presence of sampling artefacts.
}
\keyword{misc}
