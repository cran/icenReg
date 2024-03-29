% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ic_sp.R
\name{makeCtrls_icsp}
\alias{makeCtrls_icsp}
\title{Control Parameters for ic_sp}
\usage{
makeCtrls_icsp(useGA = T, maxIter = 10000, baseUpdates = 5, regStart = NULL)
}
\arguments{
\item{useGA}{Should constrained gradient ascent step be used?}

\item{maxIter}{Maximum iterations}

\item{baseUpdates}{number of baseline updates (ICM + GA) per iteration}

\item{regStart}{Initial values for regression parameters

 @description
Creates the control options for the \code{ic_sp} function. 
Defaults not intended to be changed for use in standard analyses.}
}
\description{
Control Parameters for ic_sp
}
\details{
The constrained gradient step, actived by \code{useGA = T}, 
is a step that was added to improve the convergence in a special case. 
The option to turn it off is only in place to help demonstrate it's utility. 

 \code{regStart} also for seeding of initial value of regression parameters. Intended for use in ``warm start" for bootstrap samples 
 and providing fixed regression parameters when calculating fit in qq-plots.
}
\author{
Clifford Anderson-Bergman
}
