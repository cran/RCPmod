\name{cooks.distance.regimix}
\alias{cooks.distance.regimix}
\title{Calculates leave-some-out statistics for a regimix object, principally a version of Cook's distance and cross-validated predictive logl}
\description{Performs leave-some-out measures for a regimix model. This includes a measure of how much effect leaving out an observation has on the probability of each site's RCP label. Also, this function can be used as a cross-validation workhorse.
}
\section{Method}{ cooks.distance(model, ..., oosSize = 1, times = model$n, mc.cores = 1, quiet = FALSE, type = "cooksD")
}
\arguments{
\item{ model}{a regmix object whose fit you want to assess}
\item{ ...}{ignored}
\item{ oosSize}{the size of the withheld partitions (out-of-sample size). Use 1 (default) for leave-one-out statistics, such as Cook's distance and leave-one-out validation.}
\item{ times}{the number of times to perform the re-estimation (the number of leave out groups). For each 1:times a random partition of the data, of size oosSize, is taken and the model is fitted to one of the partitions. It is predicted to the other partition. The exception is when oosSize=1 and times=model$n (leave-one-out). In such cases (the default too), the observations are left out one-by-one and not randomly.}
\item{ mc.cores}{the number of cores to spread the workload over. Default is 1. Argument is useless on Windows machines -- see \code{?parallel::mclapply}}
\item{ quiet}{should printing be suppressed? Default is no, it should not. Note that in either case, printing of the iteration trace etc is suppressed for each regimix fit.}
}
\value{
An object of class \code{regiCooksD}. It is a list of 4 elements:
\item{ Y}{the species data,}
\item{ CV}{the model$n by model$S by times array of out-of-sample predictions (this array contains a lot of NAs for where predictions would in-sample),}
\item{ cooksD}{a model$n by model$nRCP matrix of statistics that resemble Cook's distance. The statistic is the change in the prediction of RCP probability from the model with all the data to the model with only the in-sample data, and}
\item{ predLogL}{the predictive log-likelihood of each point in each withheld sample (log-likelihood contributions of withheld observations, again there will be many NAs).}
}
\seealso{\code{\link{regimix}} \code{\link{stability.regimix}}
}
\examples{
\dontrun{
#not run as R CMD check complains about the time taken.
#This code will take a little while to run (<1 minute on my computer)
#For leave-one-out cooks distance, use oosSize=1
#for serious use times will need to be larger.
system.time(\{
example( regimix);
cooksD <- cooks.distance( fm, oosSize=10, times=25)
\})
example( regimix)  #will fit a model and store in fm
#for serious use times will need to be larger.
#For leave-one-out cooks distance, use oosSize=1
cooksD <- cooks.distance( fm, oosSize=10, times=5)
}
}
\keyword{misc}
