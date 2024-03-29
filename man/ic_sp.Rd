% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ic_sp.R
\name{ic_sp}
\alias{ic_sp}
\title{Semi-Parametric models for Interval Censored Data}
\usage{
ic_sp(
  formula,
  data,
  model = "ph",
  weights = NULL,
  bs_samples = 0,
  useMCores = F,
  B = c(0, 1),
  controls = makeCtrls_icsp()
)
}
\arguments{
\item{formula}{regression formula. Response must be a \code{Surv} object of type \code{'interval2'}or \code{cbind}. See details.}

\item{data}{dataset}

\item{model}{What type of model to fit. Current choices are "\code{ph}" (Cox PH) or "\code{po}" (proportional odds)}

\item{weights}{Vector of case weights. Not standardized; see details}

\item{bs_samples}{Number of bootstrap samples used for estimation of standard errors}

\item{useMCores}{Should multiple cores be used for bootstrap sample? Does not register cluster (see example)}

\item{B}{Should intervals be open or closed? See details.}

\item{controls}{Advanced control options}
}
\description{
Fits a semi-parametric model for interval censored data. 
Can fit either a Cox-PH model or a proportional odds model.  

The covariance matrix for the regression coefficients is estimated via bootstrapping. 
For large datasets, this can become slow so parallel processing can be used to take advantage of multiple cores via the \code{foreach} package.
}
\details{
Response variable should either be of the form \code{cbind(l, u)} or 
\code{Surv(l, u, type = 'interval2')}, where \code{l} and \code{u} are the lower 
and upper ends of the interval known to contain the event of interest. 
Uncensored data can be included by setting \code{l == u}, right censored data 
can be included by setting \code{u == Inf} or \code{u == NA} and left censored 
data can be included by setting \code{l == 0}.

The argument \code{B} determines whether the intervals should be open or closed, 
i.e. \code{B = c(0,1)} implies that the event occurs within the interval \code{(l,u]}.
 The exception is that if \code{l == u}, it is assumed that the event is uncensored, 
 regardless of \code{B}. 

In regards to weights, they are not standardized. 
This means that if weight[i] = 2, this is the equivalent to having two 
observations with the same values as subject i. 

The algorithm used is inspired by the extended ICM algorithm from Wei Pan 1999.
However, it uses a conditional Newton Raphson step (for the regression parameters) 
and an ICM step (for the baseline survival parameters), rather than one single
ICM step (for both sets). In addition, a gradient ascent can also be used
to update the baseline parameters. This step is necessary if the
data contains many uncensored observations, very similar to how 
the EM algorithm greatly accelerates the ICM algorithm for the NPMLE 
(gradient ascent is used rather than the EM, as the M step is not 
in closed form for semi-parametric models). 

Earlier versions of icenReg used an active set algorithm, which was not
 as fast for large datasets.
}
\examples{
set.seed(1)

sim_data <- simIC_weib(n = 100, inspections = 5, inspectLength = 1)
ph_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, 
                data = sim_data)	
# Default fits a Cox-PH model

summary(ph_fit)		
# Regression estimates close to true 0.5 and -0.5 values


new_data <- data.frame(x1 = c(0,1), x2 = c(1, 1) )
rownames(new_data) <- c('group 1', 'group 2')
plot(ph_fit, new_data)
# plotting the estimated survival curves

po_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, 
                data = sim_data, model = 'po')
# fits a proportional odds model

summary(po_fit)

# Not run: how to set up multiple cores
# library(doParallel)
# myCluster <- makeCluster(2) 
# registerDoParallel(myCluster)
# fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2,
#              data = sim_data, useMCores = TRUE
#              bs_samples = 500)
# stopCluster(myCluster)


}
\references{
Pan, W., (1999), Extending the iterative convex minorant algorithm to the Cox model for interval-censored data, \emph{Journal of Computational and Graphical Statistics}, Vol 8(1), pp109-120

Wellner, J. A., and Zhan, Y. (1997) A hybrid algorithm for computation of the maximum likelihood estimator from censored data, \emph{Journal of the  American Statistical Association}, Vol 92, pp945-959

Anderson-Bergman, C. (preprint) Revisiting the iterative convex minorant algorithm for interval censored survival regression models
}
\author{
Clifford Anderson-Bergman
}
