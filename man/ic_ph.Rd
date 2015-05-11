\name{ic_ph}
\alias{ic_ph}
\alias{vcov.ic_coxph}
\alias{summary.ic_coxph}
\alias{plot.ic_coxph}
\title{Semi-Parametric Cox-PH model for Interval Censored Data}
\usage{
  ic_ph(formula, data, bs_samples = 20, useMCores = F, seed = NULL) 
}
\arguments{
  \item{formula}{regression formula. Response must be a \code{Surv} object of type \code{'interval2'}}
  \item{data}{dataset}
  \item{bs_samples}{Number of bootstrap samples used for estimation of standard errors}
  \item{useMCores}{Should multiple cores be used for bootstrap sample? Does not register cluster (see example)}
  \item{seed}{seed for bootstrap. If \code{seed == NULL}, a random seed is still used. See details}
}
\description{
 	Fits the semi-parametric cox-ph model for interval censored data. For the baseline survival function, a combination of ICM, VEM and nearest neighbor exchange steps are used. For regression parameters, a simple Newton-Raphson step is used. 
 	
 	The covariance matrix for the regression coefficients is estimated via bootstrapping. For large datasets, this can become slow so parallel processing can be used to take advantage of multiple cores via the \code{foreach} package. 
  }
 
\details{
	It is very important to note that a random seed is *always* set if bs_samples > 0 (via \code{set.seed(seed)}), which can create problems in simulation studies if the same seed is set in every call to ic_ph during . If \code{seed == NULL}, then the starting seed will be \code{round(runif(0, max = 10^8))}, which should be approximately equivalent to not setting a seed. 
} 
  
\examples{
	set.seed(0)

	sim_data <- simICPH_beta(n = 500)
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