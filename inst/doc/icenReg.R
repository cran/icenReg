### R code from vignette source 'icenReg.Rnw'

###################################################
### code chunk number 1: icenReg.Rnw:161-162
###################################################
library(icenReg)


###################################################
### code chunk number 2: icenReg.Rnw:165-167
###################################################
data(miceData)
head(miceData, 3)


###################################################
### code chunk number 3: icenReg.Rnw:172-173
###################################################
np_fit = ic_np(cbind(l, u) ~ grp, data = miceData)


###################################################
### code chunk number 4: icenReg.Rnw:178-180
###################################################
groupedFit1 <- ic_np(cbind(l,u) ~ 0, data = miceData)
groupedFit2 <- ic_np(miceData[,c('l', 'u')])


###################################################
### code chunk number 5: icenReg.Rnw:185-187
###################################################
plot(np_fit, col = c('blue', 'orange'),
     xlab = 'Time', ylab = 'Estimated Survival')


###################################################
### code chunk number 6: icenReg.Rnw:199-201
###################################################
data("IR_diabetes")
head(IR_diabetes, 3)


###################################################
### code chunk number 7: icenReg.Rnw:206-211
###################################################
  fit_ph <- ic_sp(cbind(left, right) ~ gender, model = 'ph', 
                  bs_samples = 100, data = IR_diabetes)
      
  fit_po <- ic_sp(cbind(left, right) ~ gender, model = 'po',
                  bs_samples = 100, data = IR_diabetes)


###################################################
### code chunk number 8: icenReg.Rnw:216-218
###################################################
  fit_po
  fit_ph


###################################################
### code chunk number 9: icenReg.Rnw:227-232
###################################################
  newdata <- data.frame(gender = c('male', 'female') )
    
  rownames(newdata) <- c('males', 'females')

  plot(fit_po, newdata)


###################################################
### code chunk number 10: icenReg.Rnw:243-245
###################################################
fit_po_gamma <- ic_par(cbind(left, right) ~ gender,
    data = IR_diabetes, model = "po", dist = "gamma")


###################################################
### code chunk number 11: icenReg.Rnw:251-252
###################################################
fit_po_gamma


###################################################
### code chunk number 12: icenReg.Rnw:257-258
###################################################
plot(fit_po_gamma, newdata, lgdLocation = "topright")


###################################################
### code chunk number 13: icenReg.Rnw:265-267
###################################################
flatPrior_fit <- ic_bayes(cbind(left, right) ~ gender,
    data = IR_diabetes, model = "po", dist = "gamma")


###################################################
### code chunk number 14: icenReg.Rnw:272-273
###################################################
flatPrior_fit


###################################################
### code chunk number 15: icenReg.Rnw:279-280
###################################################
head(flatPrior_fit$samples)  


###################################################
### code chunk number 16: icenReg.Rnw:287-288
###################################################
head(flatPrior_fit$logPosteriorDensities)  


###################################################
### code chunk number 17: icenReg.Rnw:294-296
###################################################
  plot(flatPrior_fit, newdata,
       main = 'MAP Estimates')


###################################################
### code chunk number 18: icenReg.Rnw:301-302
###################################################
  plot(flatPrior_fit$samples)


###################################################
### code chunk number 19: icenReg.Rnw:309-326
###################################################
 logPriorFunction <- function(x){
   ans <- 0 
   ans <- ans + dnorm(x[1], sd = 0.1, log = T)
   # Tight prior about 1st parameter, log_shape
   ans <- ans + dnorm(x[2], sd = 10, log = T)
   # Diffuse prior about 2nd parameter, log_scale
   ans <- ans + dnorm(x[3], sd = 0.1, log = T)
   # Tight prior about 3rd parameter, regression parameter
   return(ans)
 }    
  
informPrior_fit <- ic_bayes(cbind(left, right) ~ gender,
    data = IR_diabetes, model = "po", dist = "gamma",
    logPriorFxn = logPriorFunction)
# Fitting model with prior. 

informPrior_fit


###################################################
### code chunk number 20: icenReg.Rnw:333-341
###################################################
weak_data <- IR_diabetes[1:2,]
weakData_fit <- ic_bayes(cbind(left, right) ~ gender,
    data = weak_data,
    model = "po", dist = "gamma",
    logPriorFxn = logPriorFunction,
    controls = bayesControls(useMLE_start = F))

plot(weakData_fit$samples)


###################################################
### code chunk number 21: icenReg.Rnw:361-367
###################################################
diag_baseline(cbind(left, right) ~ gender,
    model = "po",
    data = IR_diabetes,
    dists = c("exponential", "weibull", 
              "loglogistic", "gamma"),
    lgdLocation = "topright")


###################################################
### code chunk number 22: icenReg.Rnw:372-376
###################################################
diag_baseline(fit_po, lgdLocation = "topright",
              dists = c("exponential", "weibull", 
                        "loglogistic", "gamma")
              )


###################################################
### code chunk number 23: icenReg.Rnw:395-399
###################################################
diag_covar(fit_po, lgdLocation = "topright", 
           main = "Checking Proportional Odds")
diag_covar(fit_ph, lgdLocation = "topright", 
           main = "Checking Proportional Hazards")


###################################################
### code chunk number 24: icenReg.Rnw:409-412
###################################################
diag_covar(fit_po, lgdLocation = "topright", 
           main = "Checking Proportional Odds")



###################################################
### code chunk number 25: icenReg.Rnw:414-416
###################################################
diag_covar(fit_ph, lgdLocation = "topright", 
           main = "Checking Proportional Hazards")


