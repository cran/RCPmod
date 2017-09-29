\name{plot.regimix}
\alias{plot.regimix}
\title{Plots residuals from a regimix x.}
\description{Plots the residuals from a regimix x (after calculating them).
}
\section{Method}{ plot( x, ..., type="RQR", nsim = 100, alpha.conf = c(0.9, 0.95, 0.99), quiet=FALSE)
}
\arguments{
\item{ x}{an x obtained from fitting a RCP mixture model. Such as that generated from a call to regimix(qv).}
\item{ type}{the type of residual to be plotted. Options are the default "RQR" for randomised quantile residuals (see Dunn and Smyth (1996) and Foster et al (in prep) for details) and "deviance" for the square root of minus two times the log likelihood contributions for each site (see Foster et al, 2013).}
\item{ nsim}{for type=="RQR" this argument is ignored. For type=="deviance" gives the number of simulations to use for the confidence interval. The default is 100, serious usage is likely to require more.}
\item{ alpha.conf}{the confidence level(s) to use in the residual plots for type=="deviance". Default is c(0.90,0.95,0.99). Ignored if type=="RQR."}
\item{ quiet}{should printing be performed? quiet=FALSE (default) says yes!}
\item{ species}{which species should be included in the residual plot. Default is "AllSpecies" and any subset of these (see x$names$spp for those available) is accepted.}
\item{ fitted.scale}{which scale to plot the fitted values on the x-axis? Options are "response" (default), "log" (useful sometimes for log-link models), and "logit" (useful sometimes for logit-link models).}
\item{ ...}{ignored}
}
\details{
The two types of residuals are inherently different. The "RQR" residuals produce a residual for each species at each site and the "deviance" residuals produce a site residual (no species level residual). The plots also differ, the "RQR" type generates a single normal QQ-plot for all species and all sites, and a residual versus fitted plot for all species and sites (Described in Foster et al, 2013). The "deviance" type generates a pair of Tukey mean-difference plots, similar in spirit to a QQ-plot. The first is for point-wise confidence intervals and the second is for approximate global intervals. See Foster et al (2013) for details.

The distribution for the "RQR" residuals should be standard normal. For "deviance" residuals, the distribution is unknown and simulation is used to graphically assess how odd the observed residuals look compared to ones generated assuming the model is correct.
}
\author{Scott D. Foster}
\references{
Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.

Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245

Foster, S.D., Lyons, M. and Hill, N. (in prep.) Ecological Groupings of Sample Sites in the presence of sampling artefacts.
}
\keyword{misc}
