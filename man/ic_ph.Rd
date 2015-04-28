\name{ic_ph}
\alias{ic_ph}
\alias{vcov.ic_coxph}
\alias{summary.ic_coxph}
\alias{plot.ic_coxph}
\title{Semi-Parametric Cox-PH model for Interval Censored Data}
\usage{
  ic_ph(formula, data, bs_samples = 20, useMCores = F, seed = 0) 
}
\arguments{
  \item{formula}{regression formula. Response must be a \code{Surv} object of type \code{'interval2'}}
  \item{data}{dataset}
  \item{bs_samples}{Number of bootstrap samples used for estimation of standard errors}
  \item{useMCores}{Should multiple cores be used for bootstrap sample? Does not register cluster (see example)}
  \item{seed}{seed for bootstrap}
}
\description{
 	Fits the semi-parametric cox-ph model for interval censored data. For the baseline survival function, a combination of ICM, VEM and nearest neighbor exchange steps are used. For regression parameters, a simple Newton-Raphson step is used. 
 	
 	The covariance matrix for the regression coefficients is estimated via bootstrapping. For large datasets, this can become slow so parallel processing can be used to take advantage of multiple cores via the \code{foreach} package. 
  }
\examples{
	set.seed(0)

	sim_data <- simICPH_beta(n = 1000)
	fit <- ic_ph(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data)	
	summary(fit)		
	#Regression estimates close to true 0.5 and -0.5 values


	newdata <- data.frame(x1 = c(0,1), x2 = c(1, 1) )
	rownames(newdata) <- c('group 1', 'group 2')
	plot(fit, newdata)
	#plotting the estimated survival curves

	
	# Not run: how to set up multiple cores
	# library(doParallel)
	# myCluster <- makeCluster(2, type = 'FORK') 
	# registerDoParallel(myCluster)
	# fit <- ic_ph(Surv(l, u, type = 'interval2') ~ x1 + x2,
	#              data = sim_data, useMCores = TRUE)	
}
\author{
  Clifford Anderson-Bergman
}