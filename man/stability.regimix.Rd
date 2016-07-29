\name{stability.regimix}
\alias{stability.regimix}
\title{Diagnostic checks to see if RCP groups are stable}
\description{For increasing size of hold-out samples, cooks distance and predictive log-likelihood are calculated and optionally plotted.
}
\usage{
stability.regimix(model, oosSizeRange = NULL, times = model$n,
                        mc.cores = 1, quiet = FALSE, doPlot=TRUE)
}
\arguments{
\item{ model}{a regmix model, as obtained by the function \code{regimix}. This is the model whose stability is assessed. Model must contain titbits (see ?regimix and particular attention to the argument titbits=TRUE)}
\item{ oosSizeRange}{the size of the (successive) hold-out samples. If NULL (default), then a sequence of 10 sizes, from 1 to 0.2*model$n is used. The more numbers in this range, the slower the function will run.}
\item{ times}{the number of hold-out samples to use. If times=model$n and oosSize is 1, then the sample contains each and every site. Otherwise, it is a sample of size times from the possible combinations of possible hold-out sets.}
\item{ mc.cores}{the number of cores to farm the jobs out to}
\item{ quiet}{should the progress bar be displayed (bar for each oosSizeRange)}
\item{ doPlot}{should the plots be produced? Default is that they should be.}
}
\details{The plots produced are: 1) leave-some-out Cook's distance (see \code{\link{cooks.distance.regimix}}) against holdout sample size; and 2) the predictive log-likelihood for times sites, against the holdout sample size. In both plots, the values from the original model have been added to the plot.
}
\value{\code{stability.regimix} produces a registab object. This is a list with the oosSizeRnage, disty (the mean Cook's Distance for each subset size), nRCP, n, predlogls (log-likelihood of out-of-sample sites), logl.sites (the in-sample log-likelihood for full data set).
}
\seealso{\code{\link{regimix}}, \code{\link{cooks.distance.regimix}}, \code{\link{plot.registab}}
}
\examples{
\dontrun{
#not run as R CMD check complains about the time taken.
#This code will take a little while to run (about 3.5minutes on my computer)
system.time(\{
example( regimix);
stability.regimix( fm, oosSizeRange=seq( from=1,to=fm$n\%/\%5,length=5),
      times=fm$n, mc.cores=2, doPlot=FALSE);
\})
}
}
\keyword{misc}
