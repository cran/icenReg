% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_utilities.R
\name{ic_sample}
\alias{ic_sample}
\title{Draw samples from an icenReg model}
\usage{
ic_sample(fit, newdata = NULL, sampleType = "fullSample", samples = 5)
}
\arguments{
\item{fit}{icenReg model fit}

\item{newdata}{\code{data.frame} containing covariates. If blank, will use data from model}

\item{sampleType}{type of samples See details for options}

\item{samples}{Number of samples}
}
\description{
Samples response values from an icenReg fit conditional on covariates, but not 
censoring intervals. To draw response samples conditional on covariates and 
restrained to intervals, see \code{imputeCens}.
}
\details{
Returns a matrix of samples. Each row of the matrix corresponds with a subject with the 
 covariates of the corresponding row of \code{newdata}. For each column of the matrix, 
 the same sampled parameters are used to sample response variables. 
 
 If \code{newdata} is left blank, will provide estimates for original data set. 

 There are several options for how to sample. To get random samples without accounting
 for error in the estimated parameters \code{imputeType ='fixedParSample'} takes a 
 random sample of the response variable, conditional on the response interval, 
 covariates and estimated parameters at the MLE. Alternatively, 
 \code{imputeType = 'fullSample'} first takes a random sample of the coefficients,
 (assuming asymptotic normality for the ic_par) and then takes a random sample 
 of the response variable, conditional on the response interval, 
 covariates, and the random sample of the coefficients.
}
\examples{
simdata <- simIC_weib(n = 500)

fit <- ic_par(cbind(l, u) ~ x1 + x2,
              data = simdata)

newdata = data.frame(x1 = c(0, 1), x2 = c(1,1))

sampleResponses <- ic_sample(fit, newdata = newdata, samples = 100)
}
\author{
Clifford Anderson-Bergman
}
