\name{plot.registab}
\alias{plot.registab}
\title{Diagnostic plotting to see if RCP groups are stable}
\description{For increasing size of hold-out samples, cooks distance and predictive log-likelihood are plotted.
}
\section{Method}{plot.registab(x, y, minWidth=1, ncuts=111, ylimmo, ...)
}
\arguments{
You can put normal text in \bold{Arguments}, too, like this. Remember to indent all arguments, as below.
\item{ x}{a registab object (the output of stability.regimix)}
\item{ y}{ignored.}
\item{ minWidth}{the minimum width of the density for each hold-out size. A rectangle of minWidth (at least) will be placed, and centered, at each of the hold-out sizes.}
\item{ ncuts}{number of cuts to use to describe the distributions of predicted log-likelihood. This number is a maximum value and the actual number used will depend on the data (see the breaks argument of \code{hist})}
\item{ ylimmo}{y-limits on the predicted log-likelihood plot. Typically, upper bound will be zero. If NULL (default), the lower bound is taken to be the smallest observation in the data.}
\item{ ...}{ignored}
}
\value{Nothing, but causes plots to be produced on the graphics device
}
\seealso{\code{\link{stability.regimix}}, \code{\link{regimix}}, \code{\link{cooks.distance.regimix}}
}
\examples{
\dontrun{
#not run as R CMD check complains about the time taken.
#This code will take a little while to run (about 3.5minutes on my computer)
system.time(\{
example( regimix);
my.registab <- stability.regimix( fm, oosSizeRange=seq( from=1,to=fm$n\%/\%5,length=5),
      times=fm$n, mc.cores=2, doPlot=FALSE);
plot( my.registab, minWidth=1, ncuts=15);
\})
}
}
\keyword{misc}
