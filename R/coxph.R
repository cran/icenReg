
findMaximalIntersections <- function(lower, upper){
	allVals <- sort(unique(c(lower,upper)) )
	isLeft <- allVals %in% lower
	isRight <- allVals %in% upper
	miList <- .Call("findMI", allVals, isLeft, isRight, lower, upper)
	names(miList) <- c('l_inds', 'r_inds', 'mi_l', 'mi_r')	
	return(miList)
}


fit_ICPH <- function(obsMat, covars){
	mi_info <- findMaximalIntersections(obsMat[,1], obsMat[,2])
	k = length(mi_info[['mi_l']])
	pmass <- rep(1/k, k)
	covars <- as.matrix(covars)
	myFit <- .Call('test_nne', pmass, mi_info$l_inds, mi_info$r_inds, covars) 
	names(myFit) <- c('p_hat', 'coefficients', 'final_llk', 'iterations', 'score')
	myFit[['T_bull_Intervals']] <- rbind(mi_info[['mi_l']], mi_info[['mi_r']])
	myFit$p_hat <- myFit$p_hat / sum(myFit$p_hat) 
	#removes occasional numerical error. Error on the order of 10^(-15), but causes problem when calculating last 
	#entry for estimates survival curve
	return(myFit)
}



ic_ph <- function(formula, data, bs_samples = 20, useMCores = F, seed = 0){
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    
    if (is.empty.model(mt)){
    	stop('no covariates included. Try using computeMLE in MLECens package')
    }
     x <- model.matrix(mt, mf, contrasts)
	if('(Intercept)' %in% colnames(x)){	
		ind = which(colnames(x) == '(Intercept)')
		x <- x[,-ind]
	}
	
    yMat <- as.matrix(y)[,1:2]
    rightCens <- mf[,1][,3] == 0
	yMat[rightCens,2] <- Inf
	
	exact <- mf[,1][,3] == 1
	yMat[exact, 2] = yMat[exact, 1]
    storage.mode(yMat) <- 'double'
    
    if(sum(is.na(mf)) > 0)
    	stop("NA's not allowed. If this is supposed to be right censored (i.e. [4, NA] was supposed to be right censored at t = 4), replace NA with Inf")
        
    testMat <- cbind(x, 1)
    invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
    if(is(invertResult, 'try-error'))
	    stop('covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level')
   	fitInfo <- fit_ICPH(yMat, x)
   	if(seed < 0) stop('seed must be non-negative')   
	dataEnv <- list()
	dataEnv[['x']] <- as.matrix(x, nrow = nrow(yMat))
	if(ncol(dataEnv$x) == 1) colnames(dataEnv[['x']]) <- as.character(cl[[2]][[3]])
	dataEnv[['y']] <- yMat
	seeds = 1:bs_samples + seed
	bsMat <- numeric()
	if(bs_samples > 0){
	   	if(useMCores == F){
	 		for(i in 1:bs_samples){
	    		set.seed(i + seed)
	    		sampDataEnv <- bs_sampleData(dataEnv)
				bsMat <- rbind(bsMat, getBS_coef(sampDataEnv))
	    	}
	    }
	    else{
	    	bsMat <- foreach(i = seeds, 
	    					.combine = 'rbind') %dopar%{
	    		set.seed(i)
	    		sampDataEnv <- bs_sampleData(dataEnv)
				getBS_coef(sampDataEnv)
	    	}
	    }
	    
	 xNames <- colnames(x)
	 if(!is.matrix(x)){
	 	xNames <- as.character(cl[[2]][[3]])
   	   	names(fitInfo$coefficients) <- xNames
   	 }
   	 colnames(bsMat) <- xNames
   	 incompleteIndicator <- is.na(bsMat[,1])
   	 numNA <- sum(incompleteIndicator)
   	 if(numNA > 0){
    		if(numNA / length(incompleteIndicator) >= 0.1)
    		cat('warning: ', numNA,' bootstrap samples (out of ', bs_samples, ') were dropped due to singular covariate matrix. Likely due to very sparse covariate. Be wary of these results.\n', sep = '')
    		bsMat <- bsMat[!incompleteIndicator,]
    	}
	covar <- cov(bsMat)
    est_bias <- colMeans(bsMat) - fitInfo$coefficients 
    fitInfo$coef_bc <- fitInfo$coefficients - est_bias
    }
    
    else{ 
    	bsMat <- NULL
    	covar <- NULL
    	coef_bc <- NULL
    }
    names(fitInfo$coefficients) <- xNames
    fitInfo$bsMat <- bsMat
    fitInfo$var <- covar
    fitInfo$call = cl
    fitInfo$formula = formula
    class(fitInfo) <- c('ic_coxph', 'semi_par_icph')
   return(fitInfo)
}


vcov.ic_coxph <- function(object,...) object$var

expandX <- function(formula, data){
	   cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    if ("(Intercept)" %in% colnames(x)) {
        ind = which(colnames(x) == "(Intercept)")
        x <- x[, -ind]
    }
    x <- matrix(x, nrow= nrow(mf))
    return(x)
}

getSCurves <- function(fit, newdata){
	if(inherits(fit, 'impute_par_icph'))	stop('getSCurves currently not supported for imputation model')
	if(missing(newdata)){
		grpNames <- 'baseline'
		etas = 1
	}
	else{
		if(!is.data.frame(newdata))	stop('newdata must be a data.frame with colnames = covariate names in original formula')
		grpNames <- rownames(newdata)
		reducFormula <- fit$formula
		reducFormula[[2]] <- NULL
		new_x <- expandX(reducFormula, newdata)
		log_etas <- as.numeric( new_x %*% fit$coefficients)
		etas <- exp(log_etas)
	}
	x_l <- fit$T_bull_Intervals[1,]
	x_u <- fit$T_bull_Intervals[2,]
	x_l <- c(x_l[1], x_l)
	x_u <- c(x_l[1], x_u)
	Tbull_intervals <- cbind(x_l,  x_u)
	colnames(Tbull_intervals) <- c('lower', 'upper')
	s <- 1 - c(0, cumsum(fit$p_hat))
	ans <- list(Tbull_ints = Tbull_intervals, "S_curves" = list())
	for(i in 1:length(etas)){
		eta <- etas[i]
		ans[["S_curves"]][[grpNames[i] ]] <- s^eta
	}
	return(ans)
}


plot.ic_coxph <- function(x, y, ...){
	if(inherits(x, 'impute_par_icph'))	stop('plot currently not supported for imputation model')
	
	curveInfo <- getSCurves(x, y)

	
	allx <- c(curveInfo$Tbull_ints[,1], curveInfo$Tbull_ints[,2])
	dummyx <- range(allx, finite = TRUE)
	dummyy <- c(0,1)
	
	plot(dummyx, dummyy, xlab = 'time', ylab = 'S(t)', ..., type = 'n')
	x_l <- curveInfo$Tbull_ints[,1]
	x_u <- curveInfo$Tbull_ints[,2]
	ss <- curveInfo$S_curves
	for(i in 1:length(ss)){
		lines(x_l, ss[[i]], col = i, type = 's')
		lines(x_u, ss[[i]], col = i, type = 's')
	}
	if(length(ss) > 1){
		grpNames <- names(ss)
		legend('topright', legend = grpNames, lwd = rep(1, length(grpNames) ), col = 1:length(ss))
	}
}

summary.ic_coxph <- function(object,...){
	sigfigs = 4
	fit <- object
	if(inherits(fit, 'semi_par_icph')) 	cat('\nSemi Parameteric Cox PH model for interval censored data\n')
	if(inherits(fit, 'impute_par_icph')) cat('\nMultiple Imputations Cox PH model for interval censored data\n')
	if(!is.null(fit$var)){
		colNames <- c('Estimate', 'Exp(Est)', 'Std. Error', 'z-value', 'p')
		coefs <- fit$coefficients
		output <- matrix(nrow = length(coefs), ncol = length(colNames))
		se <- sqrt(diag(fit$var))
		for(i in seq_along(coefs)){
			output[i, 1] <- coefs[i]
			output[i, 2] <- exp(coefs[i])
			output[i, 3] <- se[i]
			output[i, 4] <- coefs[i]/se[i]
			output[i, 5] <- 2*(1 - pnorm(abs(output[i,4])))
		}
		colnames(output) <- colNames
		rownames(output) <- names(coefs)
		output <- signif(output, sigfigs)
		cat("Call = \n")
		print(fit$call)
		cat('\n')
		print(output)
		if(inherits(fit, 'semi_par_icph')){
			cat('\nfinal llk = ', fit$final_llk, '\n')
			cat('Iterations = ', fit$iterations, '\n')
			cat('Bootstrap samples = ', nrow(fit$bsMat), '\n')
			if(nrow(fit$bsMat) < 100)
				cat('CAUTION: recommend more bootstrap samples for inference!\n')
		}
		
		if(inherits(fit, 'impute_par_icph')){
			cat('\nnumber of imputations = ', nrow(fit$imp_coef), '\n')
		}
	}
	else{
		colNames <- c('Estimate', 'Exp(Est)')
		coefs <- fit$coefficients
		output <- matrix(nrow = length(coefs), ncol = length(colNames))
		se <- sqrt(diag(fit$var))
		for(i in seq_along(coefs)){
			output[i, 1] <- coefs[i]
			output[i, 2] <- exp(coefs[i])
			}
		colnames(output) <- colNames
		rownames(output) <- names(coefs)
		output <- signif(output, sigfigs)
		cat("Call = \n")
		print(fit$call)
		print(output)
		cat('final llk = ', fit$final_llk, '\n')
		cat('Iterations = ', fit$iterations, '\n')	
		cat('Standard Errors not available. To get standard errors, rerun ic_ph with "bs_samples" > 0 (suggested at least 1000)')
	}
}


bs_sampleData <- function(rawDataEnv){
	n <- length(rawDataEnv[['y']][,1])
	sampEnv <- new.env()
	sampInds <- sample(1:n, n, replace = TRUE)
	sampEnv[['x']] <- rawDataEnv[['x']][sampInds,]
	sampEnv[['y']] <- rawDataEnv[['y']][sampInds,]
	return(sampEnv)
}

getBS_coef <- function(sampDataEnv){
	xMat <- cbind(sampDataEnv$x,1)
	invertResult <- try(diag(solve(t(xMat) %*% xMat )), silent = TRUE)
	if(is(invertResult, 'try-error')) {return( rep(NA, ncol(xMat) -1) ) }
	output <- fit_ICPH(sampDataEnv$y, sampDataEnv$x)$coefficients
	return(output)
}






# Multiple Imputations Regression Models



imputeCensoredData_exp <- function(l, u, impInfo, dist = 'web'){
	
	if(dist == 'exp'){
		rate <- impInfo$rates
		p_l <- pexp(l, rate)
		p_u <- pexp(u, rate)
		
		samp_q <- runif(length(l), p_l, p_u)
		samp_val <- qexp(samp_q, rate)
		is.inf <- samp_val == Inf
		samp_val[is.inf] <- u[is.inf]
		return(samp_val)
	}
	if(dist =='web'){
		rate <- impInfo$rates
		scale <- impInfo$scale
		p_l <- pweibull(l, shape = rate, scale = scale)
		p_u <- pweibull(u, shape = rate, scale = scale)
		samp_q <- runif(length(l), p_l, p_u)
		samp_val <- qweibull(samp_q, shape = rate, scale = scale)
		is.inf <- samp_val == Inf
		samp_val[is.inf] <- u[is.inf]
		return(samp_val)
	}
}


fullParamFit <- function(formula, data, param_y, dist = 'weibull'){
	data$param_y = param_y
	fit <- survreg(formula, data, dist = dist, x = TRUE)
	fit$var_chol <- chol(fit$var)
	return(fit)
}

fullParamFit_exp <- function(formula, data, param_y, dist = 'weibull'){
	data$param_y = param_y
	fit <- survreg(formula, data, dist = dist, x = TRUE)
	fit$var_chol <- chol(fit$var)
	return(fit)
}

simPars_fromFit <- function(fit, web = TRUE){
	means <- fit$coef
	if(web)
		means = c(means, log(fit$scale) )
	
	k <- length(means)
	sampPars <- means + fit$var_chol %*% rnorm(k)
	if(web){
		scale <- exp(sampPars[k])
		sampPars <- sampPars[-k]
	}
	indRates <- exp(-fit$x %*% sampPars)
	output <- list(rates = indRates)
	if(web)
		output$scale = scale
	return(output)
}

impute_ic_ph <- function(formula, data, imps = 100, eta = 10^-10, seed = 1, useMCores = FALSE){
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    
    if (is.empty.model(mt)){
    	stop('no covariates included. Try using computeMLE in MLECens package')
    }
     x <- model.matrix(mt, mf, contrasts)
	if('(Intercept)' %in% colnames(x)){	
		ind = which(colnames(x) == '(Intercept)')
		x <- x[,-ind]
	}
	if(!is(y, 'Surv'))	stop('response must be a "Surv" item with type = "interval2"')
	
    yMat <- as.matrix(y)[,1:2]
    rightCens <- mf[,1][,3] == 0
	yMat[rightCens,2] <- Inf
	exact <- mf[,1][,3] == 1
	yMat[exact, 2] = yMat[exact, 1]
	storage.mode(yMat) <- 'double'
		
	param_y <- mf[[1]]
	l_e0 <- param_y[,1] < eta
	r_e0 <- param_y[,2] < eta
	param_y[l_e0,1] <- eta
	param_y[r_e0,2] <- eta
	param_data  <- data.frame(mf[,-1])
	param_formula <- formula
	param_formula[[2]] <- as.name('param_y')
	paramFit <- fullParamFit_exp(param_formula, mf, param_y)
	
	
	fits_from_imputes <- list()
	obs_l <- yMat[,1]
	obs_r <- yMat[,2]
	use_n <- length(obs_l)
	one_vec <- rep(1, use_n)
	
	coxphFormula <- formula
	coxphFormula[[2]][[2]] <- as.name('sim_y')
	coxphFormula[[2]][[3]] <- as.name('one_vec')
	coxphFormula[[2]][[4]] <- 'right'
	
	mf_forImputes <- mf
	mf_forImputes <- data.frame(mf_forImputes)
	mf_forImputes[['one_vec']] <- one_vec
	if(!useMCores){
		set.seed(seed)
		for(i in 1:imps){
			sampPars <- simPars_fromFit(paramFit)
			sim_y <- imputeCensoredData_exp(obs_l, obs_r, sampPars)
			mf_forImputes[['sim_y']] <- sim_y		
			fits_from_imputes[[i]] <- coxph(coxphFormula, data = mf_forImputes)
		}
	}
	
	if(useMCores){
		fits_from_imputes <- foreach(i = 1:imps) %dopar%{
			set.seed(i + seed)
			sampPars <- simPars_fromFit(paramFit)
			sim_y <- imputeCensoredData_exp(obs_l, obs_r, sampPars)
			mf_forImputes[['sim_y']] <- sim_y
			return(coxph(coxphFormula, data = mf_forImputes) ) 
		}
	}
	
	k <- length(paramFit$coefficients)
	coef <- numeric()
	ave_var <- matrix(0, nrow = k-1, ncol = k-1)
	for(i in 1:imps){
		thisFit <- fits_from_imputes[[i]]
		coef <- rbind(coef, thisFit$coefficients)
		ave_var <- ave_var + thisFit$var/imps	
	}
	var <- ave_var + cov(coef)
	mean_coef <- colMeans(coef)
	fit <- list(coefficients = mean_coef, var = var, imp_coef = coef, average_imp_covariace = ave_var, call = cl, formula = formula)
    class(fit) <- c('ic_coxph', 'impute_par_icph')
	return(fit)
}

simRawTimes <- function(b1 = 0.5, b2 = -0.5, n = 100, shape1 = 2, shape2 = 2){
	rawQ <- runif(n)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    nu <- exp(x1 * b1 + x2 * b2)
    adjQ <- 1 - rawQ^(1/nu)
    trueTimes <- qbeta(adjQ, shape1 = shape1, shape2 = shape2)
	return(data.frame(y = trueTimes, x1 = x1, x2 = x2, obs = rep(1, n)))
}

simRawExpTimes <- function(b1 = 0.5, b2 = -0.5, n = 100, rate = 1){
	rawQ <- runif(n)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    nu <- exp(x1 * b1 + x2 * b2)
    adjQ <- 1 - rawQ^(1/nu)
    trueTimes <- qexp(adjQ, rate = rate)
	return(data.frame(y = trueTimes, x1 = x1, x2 = x2, obs = rep(1, n)))
}

simICPH_beta <- function(n = 100, b1 = 0.5, b2 = -0.5, inspections = 1, shape1 = 2, shape2 = 2){
	rawQ <- runif(n)
    x1 <- rnorm(n)
    x2 <- rbinom(n, 1, 0.5) - 0.5
    nu <- exp(x1 * b1 + x2 * b2)
    adjQ <- 1 - rawQ^(1/nu)
    trueTimes <- qbeta(adjQ, shape1 = shape1, shape2 = shape2)
    
    inspectionError = 1 / (inspections + 1)
    obsTimes <- 1 / (inspections + 1) + runif(n, min = -inspectionError, max = inspectionError)
    
    l <- rep(0, n)
    u <- rep(0, n)
    
    caught <- trueTimes < obsTimes
    u[caught] <- obsTimes[caught]
    l[!caught] <- obsTimes[!caught]
    
    if(inspections > 1){
    	for(i in 2:inspections){
		    oldObsTimes <- obsTimes
    		obsTimes <- i / (inspections+1) + runif(n, min = -inspectionError, max = inspectionError)
    		caught <- trueTimes >= oldObsTimes  & trueTimes < obsTimes
    		needsCatch <- trueTimes > obsTimes
    		u[caught] <- obsTimes[caught]
    		l[needsCatch] <- obsTimes[needsCatch]
    	}
    }
    else{
    	needsCatch <- !caught	
    }
    u[needsCatch] <- 1
    
    if(sum(l > u) > 0)	stop('warning: l > u! Bug in code')
    
    return(data.frame(l = l, u = u, x1 = x1, x2 = x2))
}
