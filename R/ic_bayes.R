#' Bayesian Regression  Models for Interval Censored Data
#' 
#' @param formula        Regression formula. Response must be a \code{Surv} object of type
#'  \code{'interval2'} or \code{cbind}. See details.
#' @param data           Dataset
#' @param model          What type of model to fit. Current choices are "\code{ph}" (proportional hazards), 
#' "\code{po}" (proportional odds) or "\code{aft}" (accelerated failure time)
#' @param dist           What baseline parametric distribution to use. See details for current choices
#' @param weights        vector of case weights. Not standardized; see details
#' @param logPriorFxn    An R function that computes the log prior
#' @param controls       Control parameters passed to samplers 
#'
#' @description Fits a Bayesian regression model for interval censored data. 
#' Can fita proportional hazards, proportional odds or accelerated failure time model.  
#'
#' @details Currently supported distributions choices are "exponential", "weibull", "gamma", 
#' "lnorm", "loglogistic" and "generalgamma" (i.e. generalized gamma distribution). 
#'
#' The \code{logPriorFxn} should take in the a vector of values corresponding to \emph{all}
#' the parameters of the model (baseline parameters first, regression parameters second) and returns the
#' log prior, calculated up to an additive constant. Default behavior is to use a flat prior. 
#' See examples for an example of using the log prior function.
#'
#' Sampling is done by a single MH block updater on all the parameters. 
#' See \code{?bayesControls} for more details. 
#'
#' Response variable should either be of the form \code{cbind(l, u)} or \code{Surv(l, u, type = 'interval2')}, 
#' where \code{l} and \code{u} are the lower and upper ends of the interval known to contain the event of interest. 
#' Uncensored data can be included by setting \code{l == u}, right censored data can be included by setting 
#' \code{u == Inf} or \code{u == NA} and left censored data can be included by setting \code{l == 0}.
#'
#' Does not allow uncensored data points at t = 0 (i.e. \code{l == u == 0}), as this will 
#' lead to a degenerate estimator for most parametric families. Unlike the current implementation 
#' of survival's \code{survreg}, does allow left side of intervals of positive length to 0 and 
#' right side to be \code{Inf}. 
#'
#' In regards to weights, they are not standardized. This means that if weight[i] = 2, 
#' this is the equivalent to having two observations with the same values as subject i. 
#' 
#' 
#' For numeric stability, if abs(right - left) < 10^-6, observation are considered 
#' uncensored rather than interval censored with an extremely small interval. 
#' @examples
#' data(miceData)
#' 
#' flat_prior_model <- ic_bayes(cbind(l, u) ~ grp, data = miceData)
#' # Default behavior is flat prior
#' 
#' priorFxn <- function(pars){
#'  ans <- 0
#'  ans <- ans + dnorm(pars[1], log = TRUE)
#'  ans <- ans + dnorm(pars[3], sd = 0.25, log = TRUE)
#' }
#' # Prior function puts N(0,1) prior on baseline shape parameter (first parameter)
#' # flat prior on baseline scale parameter (second parameter)
#' # and N(0,0.25) on regression parameter (third parameter)
#' 
#' inform_prior_fit <- ic_bayes(cbind(l, u) ~ grp, 
#'                              data = miceData,
#'                              logPriorFxn = priorFxn)
#' 
#' summary(flat_prior_model)
#' summary(inform_prior_fit)
#' # Note tight prior on the regression pulls posterior mean toward 0
#' 
#' @author Clifford Anderson-Bergman
#' @export
ic_bayes <- function(formula, data, logPriorFxn = function(x) return(0),
                     model = 'ph', dist = 'weibull', weights = NULL,
                     controls = bayesControls()){
  if(missing(data)) data <- environment(formula)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  #    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf, contrasts)
  if(is.matrix(x))	xNames <- colnames(x)
  else				xNames <- as.character(formula[[3]])
  if('(Intercept)' %in% colnames(x)){	
    ind = which(colnames(x) == '(Intercept)')
    x <- x[,-ind]
    xNames <- xNames[-ind]
  }
  
  yMat <- as.matrix(y)[,1:2]
  
  if(is(y, "Surv")){
    rightCens <- mf[,1][,3] == 0
    yMat[rightCens,2] <- Inf
    
    exact <- mf[,1][,3] == 1
    yMat[exact, 2] = yMat[exact, 1]
  }
  storage.mode(yMat) <- 'double'
  
  if(sum(is.na(mf)) > 0)
    stop("NA's not allowed. If this is supposed to be right censored (i.e. [4, NA] was supposed to be right censored at t = 4), replace NA with Inf")
  
  testMat <- cbind(x, 1)
  invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
  if(is(invertResult, 'try-error') & controls$useMLE_start){
    errorMsg <- 'covariate matrix is computationally singular!'
    errorMsg <- paste0(errorMsg, '\nic_bayes can still work with informative priors if useMLE_start = F in controls.')
    errorMsg <- paste0(errorMsg, '\nSee ?bayesControls for more details')
    stop( errorMsg )
  }
  modelName <- paste(dist, model, 'bayes')
  callText <- cl
  
  if(is.null(weights))	weights = rep(1, nrow(yMat))
  if(length(weights) != nrow(yMat))	stop('weights improper length!')
  if(min(weights) < 0)				stop('negative weights not allowed!')
  if(sum(is.na(weights)) > 0)			stop('cannot have weights = NA')
  if(is.null(ncol(x))) recenterCovar = FALSE
  ans <- fit_bayes(yMat, x, parFam = dist, link = model, 
                   leftCen = 0, rightCen = Inf, uncenTol = 10^-6, 
                   regnames = xNames, weights = weights,
                   callText = callText, logPriorFxn = logPriorFxn,
                   bayesList = controls, modelName = modelName)
  ans$model = model
  ans$terms <- mt
  ans$xlevels <- .getXlevels(mt, mf)
  ans$formula <- formula
  dataEnv <- new.env()
  dataEnv$data <- data
  ans$.dataEnv <- dataEnv
  return(ans)
}

#' Control parameters for ic_bayes
#' 
#' @param samples               Number of samples. 
#' @param useMLE_start          Should MLE used for starting point? 
#' @param burnIn                Number of samples discarded for burn in
#' @param samplesPerUpdate      Number of iterations between updates of proposal covariance matrix
#' @param initSD                If \code{useMLE_start == FALSE}, initial standard deviation used 
#' @param updateChol            Should cholesky decomposition be updated?
#' @param acceptRate            Target acceptance rate
#' @param thin                  Amount of thinning
#' 
#' @details 
#' 
#' Control parameters for the MH block updater used by \code{ic_bayes}.
#' 
#' The \code{samples} argument dictates how many MCMC samples are taken. One 
#' sample will be saved every \code{thin} iterations, so there will a total of
#' \code{thin * samples + burnIn} iterations. The burn in samples are not saved at all. 
#' 
#' Default behavior is to first calculate the MLE (not the MAP) estimate and use 
#' Hessian at the MLE to seed the proposal covariance matrix. After this, an updative 
#' covariance matrix is used. In cases with weakly informative likelihoods, 
#' using the MLE startpoint may lead to overly diffuse proposal or even undefined 
#' starting values. 
#' In this case, it suggested to use a cold start by setting \code{useMLE_start = F}
#' for the \code{controls} argument. In this case, the initial starting proposal
#'  covariance matrix will be a diagonal matrix with \code{initSD} standard deviations. 
#' 
#' @export
bayesControls <- function(samples = 4000, 
                          useMLE_start = TRUE, burnIn = 2999, 
                          samplesPerUpdate = 1000, initSD = 0.1,
                          updateChol = TRUE, acceptRate = 0.44,
                          thin = 5){
  ans <- list(useMLE_start        = useMLE_start,
              samples             = samples,
              thin                = thin,
              samplesPerUpdate    = samplesPerUpdate,
              updateChol          = updateChol,
              initSD              = initSD,
              burnIn              = burnIn,
              acceptRate          = acceptRate
  )
  return(ans)
}


fit_bayes <- function(y_mat, x_mat, parFam, link, 
                      leftCen = 0, rightCen = Inf, 
                      uncenTol = 10^-6, regnames, 
                      weights, callText, logPriorFxn,
                      bayesList, modelName){
  
  
  parList<- make_par_fitList(y_mat, x_mat, parFam = parFam, 
                             link = link, leftCen = 0, rightCen = Inf,
                             uncenTol = 10^-6, regnames = regnames,
                             weights = weights, callText = modelName)
  
  c_fit                  = R_ic_bayes(bayesList, logPriorFxn, parList)
  allParNames            = c(parList$bnames, parList$regnames)
  mat_samples            = c_fit$samples 
  colnames(mat_samples)  = allParNames
  mcmc_samples           = coda::mcmc(mat_samples, 
                                thin = bayesList$thin, 
                                start = bayesList$burnIn + 1)
  logPostDens            = c_fit$logPosteriorDensity
  colnames(mcmc_samples) = allParNames
  
  nBase             <- length(parList$bnames)
  nRegPar           <- length(parList$regnames)
  
  fit <- new(modelName)
  fit$par           <- parFam
  fit$baseline      <- icr_colMeans(mat_samples[ ,1:nBase])
  regParVec         <- NULL
  covMat            <- NULL
  if(nRegPar > 0) { 
    regParVec <-  icr_colMeans(mat_samples[ ,nBase + 1:nRegPar])
  } 
  fit$reg_pars      <- regParVec
  fit$nSamples      <- nrow(mat_samples)
  fit$var           <- cov(mat_samples)
  fit$samples       <- mcmc_samples
  fit$logPosteriorDensities <- logPostDens
  fit$ess           <- coda::effectiveSize(mcmc_samples)
  fit$call          <- callText
  fit$logPrior      <- logPriorFxn
  fit$finalChol     <- c_fit$finalChol
  fit$MAP_ind       <- which(fit$logPosteriorDensities == max(fit$logPosteriorDensities) )[1]
  fit$MAP_reg_pars  <- fit$samples[fit$MAP_ind, nBase + 1:nRegPar]
  fit$MAP_baseline  <- fit$samples[fit$MAP_ind, 1:nBase]
  names(fit$reg_pars)    <- parList$regnames
  names(fit$baseline)    <- parList$bnames
  names(fit$MAP_reg_pars)    <- parList$regnames
  names(fit$MAP_baseline)    <- parList$bnames
  fit$coefficients <- c(fit$baseline, fit$reg_pars)
  fit$baseOffset   <- 0
  return(fit)
}