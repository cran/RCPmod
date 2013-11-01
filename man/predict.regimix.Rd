\name{predict.regimix}
\alias{predict.regimix}
\title{Predicts from a regimix object.}
\description{Predicts RCP probabilities at a series of sites. Confidence intervals are available too.
}
\section{Method}{ predict( object, ..., newdata=NULL, nboot=0, alpha=0.95)
}
\arguments{
\item{ object}{an object obtained from fitting a RCP mixture model. Such as that generated from a call to regimix(qv).}
\item{ newdata}{a data.frame (or something that can be coerced) containing the values of the covariates where predictions are to be made. If NULL (the default) then predictions are made at the locations of the original data.}
\item{ nboot}{the number of parametric bootstrap samples to take for the bootstrap predictions, standard errors and confidence intervals. The default is 0, that is no bootstrapping is to be done and point predictions only are given.}
\item{ alpha}{a numeric within [0,1] (well [0.5,1] really) indicating the specified confidence for the confidence interval. Argument is redundant if nboot == 0.}
\item{ ...}{ignored}
}
\details{This function can take a while to run -- it is a bootstrap function. nboot samples of the parameters are taken and then used to predict at each set of covariates defined in newdata. Quantiles of the resulting sets of bootstrap predictions are then taken. It is the last step that really takes a while.
}
\value{
If nboot==0 then a n x H matrix of prior predictions (n=nrow(newdata), H=number of RCPs). Each row should sum to one.
\item{ if nboot!=0 then a list is returned. It has elements:}{}
\item{ ptPreds}{the n x H matrix of point predictions}
\item{ bootPreds}{the n x H matrix of bootstrap point predictions (mean of bootstrap samples)}
\item{ bootSEs}{the n x H matrix of bootstrap standard errors for predictions}
\item{ bootCIs}{the n x H x 2 array of bootstrap confidence intervals. Note that bootCIs[,,1] gives the lower CIs and bootCIs[,,2] gives the upper CIs.}
}
\author{Scott D. Foster}
\keyword{misc}
