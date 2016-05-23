### R code from vignette source 'icenReg.Rnw'

###################################################
### code chunk number 1: icenReg.Rnw:155-156
###################################################
library(icenReg)


###################################################
### code chunk number 2: icenReg.Rnw:159-161
###################################################
data(miceData)
head(miceData, 3)


###################################################
### code chunk number 3: icenReg.Rnw:166-167
###################################################
np_fit = ic_np(cbind(l, u) ~ grp, data = miceData)


###################################################
### code chunk number 4: icenReg.Rnw:172-174
###################################################
groupedFit1 <- ic_np(cbind(l,u) ~ 0, data = miceData)
groupedFit2 <- ic_np(miceData[,c('l', 'u')])


###################################################
### code chunk number 5: icenReg.Rnw:179-181
###################################################
plot(np_fit, col = c('blue', 'orange'),
     xlab = 'Time', ylab = 'Estimated Survival')


###################################################
### code chunk number 6: icenReg.Rnw:193-195
###################################################
data("IR_diabetes")
head(IR_diabetes, 3)


###################################################
### code chunk number 7: icenReg.Rnw:200-205
###################################################
  fit_ph <- ic_sp(cbind(left, right) ~ gender, model = 'ph', 
                  bs_samples = 100, data = IR_diabetes)
      
  fit_po <- ic_sp(cbind(left, right) ~ gender, model = 'po',
                  bs_samples = 100, data = IR_diabetes)


###################################################
### code chunk number 8: icenReg.Rnw:210-212
###################################################
  fit_po
  fit_ph


###################################################
### code chunk number 9: icenReg.Rnw:221-226
###################################################
  newdata <- data.frame(gender = c('male', 'female') )
    
  rownames(newdata) <- c('males', 'females')

  plot(fit_po, newdata)


###################################################
### code chunk number 10: icenReg.Rnw:237-239
###################################################
fit_po_gamma <- ic_par(cbind(left, right) ~ gender,
    data = IR_diabetes, model = "po", dist = "gamma")


###################################################
### code chunk number 11: icenReg.Rnw:245-246
###################################################
fit_po_gamma


###################################################
### code chunk number 12: icenReg.Rnw:251-252
###################################################
plot(fit_po_gamma, newdata, lgdLocation = "topright")


###################################################
### code chunk number 13: icenReg.Rnw:270-276
###################################################
diag_baseline(cbind(left, right) ~ gender,
    model = "po",
    data = IR_diabetes,
    dists = c("exponential", "weibull", 
              "loglogistic", "gamma"),
    lgdLocation = "topright")


###################################################
### code chunk number 14: icenReg.Rnw:281-285
###################################################
diag_baseline(fit_po, lgdLocation = "topright",
              dists = c("exponential", "weibull", 
                        "loglogistic", "gamma")
              )


###################################################
### code chunk number 15: icenReg.Rnw:304-308
###################################################
diag_covar(fit_po, lgdLocation = "topright", 
           main = "Checking Proportional Odds")
diag_covar(fit_ph, lgdLocation = "topright", 
           main = "Checking Proportional Hazards")


###################################################
### code chunk number 16: icenReg.Rnw:318-321
###################################################
diag_covar(fit_po, lgdLocation = "topright", 
           main = "Checking Proportional Odds")



###################################################
### code chunk number 17: icenReg.Rnw:323-325
###################################################
diag_covar(fit_ph, lgdLocation = "topright", 
           main = "Checking Proportional Hazards")


