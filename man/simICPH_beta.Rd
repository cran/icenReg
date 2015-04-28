\name{simICPH_beta}
\alias{simICPH_beta}
\title{Simulates interval censored data from Cox-PH model with beta baseline}
\usage{
  simICPH_beta(n = 100, b1 = 0.5, b2 = -0.5, inspections = 1, shape1 = 2, shape2 = 2)
}
\arguments{
  \item{n}{Number of samples simulated}
  \item{b1}{Value of first regression coefficient}
  \item{b2}{Value of second regression coefficient}
  \item{inspections}{Number of inspections of data. inspections = 1 implies current status data}
  \item{shape1}{Shape1 parameter of baseline beta distribution}
  \item{shape2}{Shape2 parameter of baseline beta distribution}
}
\description{
  
  Simulates interval censored data from a Cox-PH model with a beta baseline distribution. Used for demonstration.}
\examples{
	set.seed(0)

	sim_data <- simICPH_beta(n = 1000)
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

	sim_data_fine <- simICPH_beta(n = 1000, inspections = 20)	
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
	#will be negative the same coefficients for a cox-ph)}
\author{
  Clifford Anderson-Bergman
}