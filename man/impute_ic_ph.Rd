\name{impute_ic_ph}
\alias{impute_ic_ph}
\title{Multiple Imputations Cox-PH model for Interval Censored Data}
\usage{
  impute_ic_ph(formula, data, imps = 100, 
               eta = 10^-10, rightCenVal = 10000,
               seed = NULL, useMCores = FALSE) 
}
\arguments{
  \item{formula}{regression formula. Response must be a \code{Surv} object of type \code{'interval2'}}
  \item{data}{dataset}
  \item{imps}{Number of imputations used}
  \item{eta}{0 offset. See details}
  \item{rightCenVal}{Value to be used for right side of right censored data}
  \item{seed}{seed for random number generator. If \code{seed == NULL}, a random seed is still used. See details}
  \item{useMCores}{Should multiple cores be used for imputation? See \code{?ic_ph} for example of how to set up}
}
\description{
	Analyzes interval censored data via multiple imputation. In particular, it first fits a fully parametric model with an Weibull baseline distribution (via \code{survival}'s \code{survreg} function). Then, it uses this model to impute the exact times for each censored interval. These imputed datasets are then fit with a semi-parametric model (via \code{survival}'s \code{coxph} function). 
	
	In general, it is much better to use the \code{ic_ph} function. However, for very large datasets, this may be considerably faster. And if the length of the intervals is relatively small, then the estimators should agree heavily with each other. However, for data with larger intervals, such as current status data, the assumption of a Weibull baseline distribution can be influential. See examples for demonstration. 
}
\details{The survival function \code{survreg}, which impute_ic_ph uses for the imputation model, does not allow for intervals to start with 0, nor intervals that end with Inf. \code{impute_ic_ph} adds \code{eta} to left censored data to prevent this from being an issue and replaces the right hand side of an interval that ends with \code{Inf} with \code{rightCenVal}.
	
		It is very important to note that a random seed is *always* set (via \code{set.seed(seed)}), which can create problems in simulation studies if the same seed is set in every call to impute_ic_ph during . If \code{seed == NULL}, then the starting seed will be \code{round(runif(0, max = 10^8))}, which should be approximately equivalent to not setting a seed. 
}
\examples{
	set.seed(0)

	sim_data <- simICPH_beta(n = 500)
	#simulates data from a cox-ph with beta baseline distribution.
	#default inspections = 1 results in current status data
	
	
	fit <- ic_ph(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data)	
	summary(fit)		
	#Analyzing with semi-parametric model
	#Regression estimates close to true 0.5 and -0.5 values
	
	imp_fit <- impute_ic_ph(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data)
	summary(imp_fit)		
	#Analyzing with imputation cox ph model
	#Heavy bias from grossly inapporiate assumption of baseline for imputations 
	#(true baseline = beta, baseline used for imputations = Weibull)

	sim_data_fine <- simICPH_beta(n = 500, inspections = 20)	
	#Now there are more inspections
	#The fact that imputed values are restrained to narrow intervals downplays the 
	#negative influence of the Weibull imputation assumption

	imp_fit2 <- impute_ic_ph(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data_fine)
	summary(imp_fit2)
	#Results in a less biased estimator
	
	fullExpFit <- survreg(Surv(l + 10^-10, u, type = 'interval2') ~ x1 + x2,
	                      data = sim_data_fine, dist = 'exponential')
	summary(fullExpFit)
	#Note that using a fully parametric model still results in heavy bias despite fine 
	#grid of inspection times
	#(Note: survreg fits aft model, but for exponential distribution, coefficients
	#will be -1 * coefficients for a cox-ph)
}
\author{
  Clifford Anderson-Bergman
}