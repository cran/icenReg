### R code from vignette source 'icenReg.Rnw'

###################################################
### code chunk number 1: icenReg.Rnw:158-159
###################################################
library(icenReg)


###################################################
### code chunk number 2: icenReg.Rnw:162-164
###################################################
data(miceData)
head(miceData, 3)


###################################################
### code chunk number 3: icenReg.Rnw:169-170
###################################################
np_fit = ic_np(cbind(l, u) ~ grp, data = miceData)


###################################################
### code chunk number 4: icenReg.Rnw:175-177
###################################################
groupedFit1 <- ic_np(cbind(l,u) ~ 0, data = miceData)
groupedFit2 <- ic_np(miceData[,c('l', 'u')])


###################################################
### code chunk number 5: icenReg.Rnw:182-184
###################################################
plot(np_fit, col = c('blue', 'orange'),
     xlab = 'Time', ylab = 'Estimated Survival')


###################################################
### code chunk number 6: icenReg.Rnw:196-198
###################################################
data("IR_diabetes")
head(IR_diabetes, 3)


###################################################
### code chunk number 7: icenReg.Rnw:203-208
###################################################
  fit_ph <- ic_sp(cbind(left, right) ~ gender, model = 'ph', 
                  bs_samples = 100, data = IR_diabetes)
      
  fit_po <- ic_sp(cbind(left, right) ~ gender, model = 'po',
                  bs_samples = 100, data = IR_diabetes)


###################################################
### code chunk number 8: icenReg.Rnw:213-215
###################################################
  fit_po
  fit_ph


###################################################
### code chunk number 9: icenReg.Rnw:224-229
###################################################
  newdata <- data.frame(gender = c('male', 'female') )
    
  rownames(newdata) <- c('males', 'females')

  plot(fit_po, newdata)


###################################################
### code chunk number 10: icenReg.Rnw:240-242
###################################################
fit_po_gamma <- ic_par(cbind(left, right) ~ gender,
    data = IR_diabetes, model = "po", dist = "gamma")


###################################################
### code chunk number 11: icenReg.Rnw:248-249
###################################################
fit_po_gamma


###################################################
### code chunk number 12: icenReg.Rnw:254-255
###################################################
plot(fit_po_gamma, newdata, lgdLocation = "topright")


###################################################
### code chunk number 13: icenReg.Rnw:273-279
###################################################
diag_baseline(cbind(left, right) ~ gender,
    model = "po",
    data = IR_diabetes,
    dists = c("exponential", "weibull", 
              "loglogistic", "gamma"),
    lgdLocation = "topright")


###################################################
### code chunk number 14: icenReg.Rnw:284-288
###################################################
diag_baseline(fit_po, lgdLocation = "topright",
              dists = c("exponential", "weibull", 
                        "loglogistic", "gamma")
              )


###################################################
### code chunk number 15: icenReg.Rnw:307-311
###################################################
diag_covar(fit_po, lgdLocation = "topright", 
           main = "Checking Proportional Odds")
diag_covar(fit_ph, lgdLocation = "topright", 
           main = "Checking Proportional Hazards")


###################################################
### code chunk number 16: icenReg.Rnw:321-324
###################################################
diag_covar(fit_po, lgdLocation = "topright", 
           main = "Checking Proportional Odds")



###################################################
### code chunk number 17: icenReg.Rnw:326-328
###################################################
diag_covar(fit_ph, lgdLocation = "topright", 
           main = "Checking Proportional Hazards")


