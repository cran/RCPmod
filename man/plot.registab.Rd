\name{plot.registab}
\alias{plot.registab}
\title{Diagnostic plotting to see if RCP groups are stable}
\description{For increasing size of hold-out samples, cooks distance and predictive log-likelihood are plotted.
}
\section{Method}{plot.registab(x, y, ...)
}
\arguments{
You can put normal text in \bold{Arguments}, too, like this. Remember to indent all arguments, as below.
\item{ x}{a registab object (the output of stability.regimix)}
\item{ y}{ignored.}
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
plot( my.registab);
\})
}
}
\keyword{misc}
