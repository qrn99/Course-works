setwd("~/Desktop/UofT/19_20/STA302")
# install.packages('NHANES'); install.packages('Matrix'); install.packages('tidyverse'); install.packages('glmnet'); install.packages('rms')
# install.packages("funModeling"); install.packages("Hmisc"); install.packages('DataExplorer'); install.packages("broom"); install.packages("purrr"); install.packages("jtools")

## Final Project STA302 Summer 2020 ##
rm(list = ls())
library(NHANES); library(tidyverse); library(glmnet); library(rms); library(xtable); library(car); 
library(funModeling); library(Hmisc); library(DataExplorer); library(broom); library(purrr); library(jtools)  # data description and t-test, final model result output

small.nhanes <- na.omit(NHANES[NHANES$SurveyYr=="2011_12"
                               & NHANES$Age > 17,c(1,3,4,8:11,13,17,20,21,25,46,50,51,52,61)])
small.nhanes <- as.data.frame(small.nhanes %>%
  group_by(ID) %>% filter(row_number()==1) )
nrow(small.nhanes)

## Model selection criteria ##
criteria <- function(model){
  n <- length(model$residuals)
  p <- length(model$coefficients) - 1
  RSS <- sum(model$residuals^2)
  R2 <- summary(model)$r.squared
  R2.adj <- summary(model)$adj.r.squared
  AIC <- n*log(RSS/n) + 2*p
  AICc <- AIC + (2*(p+2)*(p+3))/(n-p-1)
  BIC <- n*log(RSS/n) + (p+2)*log(n)
  res <- c(R2, R2.adj, AIC, AICc, BIC)
  names(res) <- c("R Squared", "Adjsuted R Squared", "AIC", "AICc", "BIC")
  return(res)
}

## Checking whether there are any ID that was repeated. If not ##
## then length(unique(small.nhanes$ID)) and nrow(small.nhanes) are same ##
length(unique(small.nhanes$ID))

## Create training and test set ##
set.seed(1004079631)
train <- small.nhanes[sample(seq_len(nrow(small.nhanes)), size = 400), ]  # Training set
nrow(train)
length(which(small.nhanes$ID %in% train$ID))
test <- small.nhanes[!small.nhanes$ID %in% train$ID,]  # Testing set
nrow(test)

# Data description:
# data_des_0 <- summary(train[, -c(1,12)])
data_des_1 <- summary(train[,-c(1, 12, 7:17)])
data_des_2 <- summary(train[, -c(1, 12, 2:6, 12:17)])
data_des_3 <- summary(train[, -c(1, 12, 2:12)])
data_des_4 <- summary(train[, -c(2:11, 13:17)])
xtable(data_des_1)
xtable(data_des_2)
xtable(data_des_3)
xtable(data_des_4)


freq(train)
xtable(freq(train))
plot_num(train)
plot_bar(train)


## Running the MLR model ##
### First fit a multiple linear regression using the training set##
model.lm <- lm(BPSysAve ~ ., data = train[, -c(1)])
summary(model.lm); # Check p-value, R^2 adjusted
# criteria
crit_lm <- criteria(model.lm)

# Normality check
r_lm <- rstudent(model.lm)
Y_hat_lm <- predict(model.lm)
qqnorm(r_lm); qqline(r_lm)

# Standardized Residual plot for fitted y-values on all 15 predictors, Y_hat_lm
plot(Y_hat_lm, r_lm, type = "p", xlab = "Fitted Values Y_hat_lm", ylab = "Standardized Residuals", main = "Standardized Residuals Plot with Fitted Values Y_hat_lm", col = "blue")
abline(h = 2, lty = 2); abline(h = -2, lty = 2); abline(h = 0, lwd = 0.8); lines(lowess(Y_hat_lm, r_lm), col="red")

# Response vs fitted values Y_hat_lm
plot(train$BPSysAve ~ Y_hat_lm, type="p", xlab="Fitted Values Y_hat_lm", ylab="Y", main = "Scatter Plot with Y vs Fitted Values Y_hat_lm", cex.lab=1.2, col='red')
abline(lm(train$BPSysAve ~ Y_hat_lm), lwd=2, col="blue"); lines(lowess(Y_hat_lm, train$BPSysAve), col="green")

# VIF
vif_lm <- vif(model.lm); which(vif_lm > 5)

# inf_obs
# hat values
h_lm <- hatvalues(model.lm); thresh <- 2 * 36/400; which(h_lm > thresh)

# Cook's distance
D_lm <- cooks.distance(model.lm); which(D_lm > qf(0.5, 36, 400-35-1))  # 50th percentile of F with p+1 and n-p-1 df

# DFFITS
dfits_lm <- dffits(model.lm); which(abs(dfits_lm) > 2*sqrt(36/400))  # 2*sqrt((p+1)/n)

# DFBETAS
dfb_lm <- dfbetas(model.lm); which(abs(dfb_lm[,1]) > 2/sqrt(400))  # 2*sqrt(1/n)

## Perform Prediction on the test set ##
pred.y <- predict(model.lm, newdata = test, type = "response")
## Prediction error ##
mean((test$BPSysAve - pred.y)^2)

# T-test for smokenow
t.test(train$BPSysAve ~ train$SmokeNow)

## Model Selections ##
#### Based on AIC ####
sel.var.aic <- step(model.lm, trace = 0, k = 2, direction = "both") 
sel.var.aic <- attr(terms(sel.var.aic), "term.labels")   
sel.var.aic

## Diagnostics ##
model.aic <- lm(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.aic, "BPSysAve"))])
crit_aic <- criteria(model.aic)

# VIF
vif_aic <- vif(model.aic); which(vif_aic > 5)

# Normality check
r_aic <- rstudent(model.aic)
Y_hat_aic <- predict(model.aic)
qqnorm(r_aic); qqline(r_aic)

# Standardized Residual plot for fitted y-values on AIC selected varibles, Y_hat_aic
plot(Y_hat_aic, r_aic, type = "p", xlab = "Fitted Values Y_hat_aic", ylab = "Standardized Residuals", main = "Standardized Residuals Plot with Fitted Values Y_hat_aic", col = "blue")
abline(h = 2, lty = 2); abline(h = -2, lty = 2); abline(h = 0, lwd = 0.8); lines(lowess(Y_hat_aic, r_aic), col="red")

# Response vs fitted values Y_hat_aic
plot(train$BPSysAve ~ Y_hat_aic, type="p", xlab="Fitted Values Y_hat_aic", ylab="Y", main = "Scatter Plot with Y vs Fitted Values Y_hat_aic", cex.lab=1.2, col='red')
abline(lm(train$BPSysAve ~ Y_hat_aic), lwd=2, col="blue"); lines(lowess(Y_hat_aic, train$BPSysAve), col="green")

### Cross Validation and prediction performance of AIC based selection ###
ols.aic <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.aic, "BPSysAve"))], 
               x=T, y=T, model = T)

## 10 fold cross validation ##    
aic.cross <- calibrate(ols.aic, method = "crossvalidation", B = 10)
## Calibration plot ##
# pdf("aic_cross.pdf", height = 8, width = 16)
plot(aic.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with AIC Model")
dev.off()

## Test Error ##
pred.aic <- predict(ols.aic, newdata = test[,which(colnames(train) %in% c(sel.var.aic, "BPSysAve"))])
## Prediction error ##
pred.error.AIC <- mean((test$BPSysAve - pred.aic)^2)
pred.error.AIC

#### Based on BIC ####
n <- nrow(train)
sel.var.bic <- step(model.lm, trace = 0, k = log(n), direction = "both") 
sel.var.bic <- attr(terms(sel.var.bic), "term.labels")   
sel.var.bic

## Diagnostics ##
model.bic <- lm(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.bic, "BPSysAve"))])
crit_bic <- criteria(model.bic)

# VIF
vif_bic <- vif(model.bic); which(vif_bic > 5)

# Normality check
r_bic <- rstudent(model.bic)
Y_hat_bic <- predict(model.bic)
qqnorm(r_bic); qqline(r_bic)

# Standardized Residual plot for fitted y-values on BIC selected varibles, Y_hat_bic
plot(Y_hat_bic, r_bic, type = "p", xlab = "Fitted Values Y_hat_bic", ylab = "Standardized Residuals", main = "Standardized Residuals Plot with Fitted Values Y_hat_bic", col = "blue")
abline(h = 2, lty = 2); abline(h = -2, lty = 2); abline(h = 0, lwd = 0.8); lines(lowess(Y_hat_bic, r_bic), col="red")

# Response vs fitted values Y_hat_bic
plot(train$BPSysAve ~ Y_hat_bic, type="p", xlab="Fitted Values Y_hat_bic", ylab="Y", main = "Scatter Plot with Y vs Fitted Values Y_hat_bic", cex.lab=1.2, col='red')
abline(lm(train$BPSysAve ~ Y_hat_bic), lwd=2, col="blue"); lines(lowess(Y_hat_bic, train$BPSysAve), col="green")

### Cross Validation and prediction performance of BIC based selection ###
ols.bic <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.bic, "BPSysAve"))], 
               x=T, y=T, model = T)

## 10 fold cross validation ##   
set.seed(1004079631)
bic.cross <- calibrate(ols.bic, method = "crossvalidation", B = 10)
## Calibration plot ##
# pdf("bic_cross.pdf", height = 8, width = 16)
plot(bic.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with BIC Model")
dev.off()

## Test Error ##
pred.bic <- predict(ols.bic, newdata = test[,which(colnames(train) %in% c(sel.var.bic, "BPSysAve"))])
## Prediction error ##
pred.error.BIC <- mean((test$BPSysAve - pred.bic)^2)
pred.error.BIC


#### Elastic Net ####
model.EN <- glmnet(x = model.matrix( ~ ., data = train[,-c(1,12)]), y = train$BPSysAve, standardize = T, alpha = 0.5)

## Perform cross validation to choose lambda for Elastic Net##
set.seed(1004079631)
cv.out_EN <- cv.glmnet(x = model.matrix( ~ ., data = train[,-c(1,12)]), y = train$BPSysAve, standardize = T, alpha = 0.5)
plot(cv.out_EN)
best.lambda_EN <- cv.out_EN$lambda.1se
best.lambda_EN
co_EN <- coef(cv.out_EN, s = "lambda.1se")

#Selection of the significant features(predictors)
## threshold for variable selection ##
thresh <- 0.00
# select variables #
inds_EN <- which(abs(co_EN) > thresh)
variables_EN <- row.names(co_EN)[inds_EN]
sel.var.EN <- variables_EN[!(variables_EN %in% '(Intercept)')]
sel.var.EN

### Cross Validation and prediction performance of AIC based selection ###
ols.EN <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.EN, "BPSysAve"))], 
               x=T, y=T, model = T)

## 10 fold cross validation ##    
set.seed(1004079631)
EN.cross <- calibrate(ols.EN, method = "crossvalidation", B = 10)
## Calibration plot ##
# pdf("aic_cross.pdf", height = 8, width = 16)
plot(EN.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with Elastic Net Model")
dev.off()

model.EN_cv <- lm(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.EN, "BPSysAve"))])

## Diagnostics ##
crit_EN <- criteria(model.EN_cv)

# Normality check
r_EN <- rstudent(model.EN_cv)
Y_hat_EN <- predict(model.EN_cv)
qqnorm(r_EN); qqline(r_EN)

# Standardized Residual plot for fitted y-values on EN selected varibles, Y_hat_EN
plot(Y_hat_EN, r_EN, type = "p", xlab = "Fitted Values Y_hat_EN", ylab = "Standardized Residuals", main = "Standardized Residuals Plot with Fitted Values Y_hat_EN", col = "blue")
abline(h = 2, lty = 2); abline(h = -2, lty = 2); abline(h = 0, lwd = 0.8); lines(lowess(Y_hat_EN, r_EN), col="red")

# Response vs fitted values Y_hat_EN
plot(train$BPSysAve ~ Y_hat_EN, type="p", xlab="Fitted Values Y_hat_EN", ylab="Y", main = "Scatter Plot with Y vs Fitted Values Y_hat_EN", cex.lab=1.2, col='red')
abline(lm(train$BPSysAve ~ Y_hat_EN), lwd=2, col="blue"); lines(lowess(Y_hat_EN, train$BPSysAve), col="green")

## Perform Prediction ##
pred.y.EN <- predict(model.EN_cv, newdata = test[,which(colnames(train) %in% c(sel.var.EN, "BPSysAve"))], type = "response")
## Prediction error after CV for EN##
pred.error.EN <- mean((test$BPSysAve - pred.y.EN)^2)
pred.error.EN


#### LASSO ####
model.lasso <- glmnet(x = model.matrix( ~ ., data = train[,-c(1,12)]), y = train$BPSysAve, standardize = T, alpha = 1)

## Perform cross validation to choose lambda ##
set.seed(1004079631)
cv.out <- cv.glmnet(x = model.matrix( ~ ., data = train[,-c(1,12)]), y = train$BPSysAve, standardize = T, alpha = 1)
plot(cv.out)
best.lambda <- cv.out$lambda.1se
best.lambda
co <- coef(cv.out, s = "lambda.1se")

#Selection of the significant features(predictors)
## threshold for variable selection ##
thresh <- 0.00
# select variables #
inds <- which(abs(co) > thresh )
variables <- row.names(co)[inds]
sel.var.lasso <- variables[!(variables %in% '(Intercept)')]
sel.var.lasso

model.lasso_cv <- lm(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.lasso, "BPSysAve"))])

## Diagnostics ##
crit_lasso_cv <- criteria(model.lasso_cv)
crit_lasso_cv

# Normality check
r_lasso_cv <- rstudent(model.lasso_cv)
Y_hat_lasso_cv <- predict(model.lasso_cv)
qqnorm(r_lasso_cv); qqline(r_lasso_cv)

# Standardized Residual plot for fitted y-values on lasso_cv  selected varibles, Y_hat_lasso_cv 
plot(Y_hat_lasso_cv, r_lasso_cv, type = "p", xlab = "Fitted Values Y_hat_lasso_cv", ylab = "Standardized Residuals", main = "Standardized Residuals Plot with Fitted Values Y_hat_lasso_cv", col = "blue")
abline(h = 2, lty = 2); abline(h = -2, lty = 2); abline(h = 0, lwd = 0.8); lines(lowess(Y_hat_lasso_cv, r_lasso_cv), col="red")

# Response vs fitted values Y_hat_lasso_cv
plot(train$BPSysAve ~ Y_hat_lasso_cv, type="p", xlab="Fitted Values Y_hat_lasso_cv", ylab="Y", main = "Scatter Plot with Y vs Fitted Values Y_hat_lasso_cv", cex.lab=1.2, col='red')
abline(lm(train$BPSysAve ~ Y_hat_lasso_cv), lwd=2, col="blue"); lines(lowess(Y_hat_lasso_cv, train$BPSysAve), col="green")


### Cross Validation and prediction performance of lasso based selection ###
ols.lasso <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.lasso, "BPSysAve"))], 
                 x=T, y=T, model = T)

## 10 fold cross validation ##    
set.seed(1004079631)
lasso.cross <- calibrate(ols.lasso, method = "crossvalidation", B = 10)
## Calibration plot ##
# pdf("lasso_cross.pdf", height = 8, width = 16)
plot(lasso.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with LASSO")
dev.off()

## Perform Prediction for ols function##
pred.lasso <- predict(ols.lasso, newdata = test[,which(colnames(train) %in% c(sel.var.lasso, "BPSysAve"))])
## Prediction error ##
pred.error.lasso <- mean((test$BPSysAve - pred.lasso)^2)
pred.error.lasso

## Perform Prediction for lm function##
pred.y.lasso_cv <- predict(model.lasso_cv, newdata = test[,which(colnames(train) %in% c(sel.var.lasso, "BPSysAve"))], type = "response")
## Prediction error ##
pred.error.lasso_cv <- mean((test$BPSysAve - pred.y.lasso_cv)^2)
pred.error.lasso_cv  # ols and lm same result

#### summary of criteria ####
crit_aic[6] = pred.error.AIC
crit_bic[6] = pred.error.BIC
crit_EN[6] = pred.error.EN
crit_lasso_cv[6] = pred.error.lasso_cv

crit_all <- matrix(c(crit_aic ,crit_bic, crit_EN, crit_lasso_cv),ncol=6,byrow=TRUE)
colnames(crit_all) <- c("R Squared", "Adjsuted R Squared", "AIC", "AICc", "BIC", "Predcition Error")
rownames(crit_all) <- c("AIC Model", "BIC Model", "Elastic-Net Model", "LASSO Model")
crit_all <- as.table(crit_all)
crit_all
xtable(crit_all)

# Calibration plot in 2by2 margin of the 4 candidate models
pdf("calib_plot.pdf", height = 8, width = 16)
par(mfrow = c(2,2), mai = c(1, 0.5, 0.3, 0.3))
plot(aic.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with AIC Model")
plot(bic.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with BIC Model")
plot(EN.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with Elastic-Net Model")
plot(lasso.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with LASSO Model")
dev.off()

# T-test on SmokeNow
t_test <- t.test(train$BPSysAve ~ train$SmokeNow)
t_test
tab <- map_df(list(t_test), tidy)
tab_t_test <- tab[c("estimate", "statistic", "p.value", "conf.low", "conf.high")]
tab_t_test
xtable(tab_t_test)


#### Final Model ####
# BIC with smokenow
model.bic_smokenow <- lm(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.bic,"SmokeNow", "BPSysAve"))])
## Diagnostics ##
crit_bic_smokenow <- criteria(model.bic_smokenow)
crit_bic_smokenow; crit_bic
summary(model.bic_smokenow)
sum_bic_somkenow <- summ(model.bic_smokenow, digits = 4)
sum_bic_somkenow

# VIF
vif_bic_smokenow  <- vif(model.bic_smokenow); which(vif_bic_smokenow > 5)

# inf_obs
# hat values
length(model.bic_smokenow$coefficients)-1 # the number of parameters that we fitted excluding the intercept as it is a separate parameter, the value of p
h_bic_smokenow <- hatvalues(model.bic_smokenow); thresh <- 2 * 13/400; which(h_bic_smokenow > thresh)  # 2*(p+1)/n

# Cook's distance
D_bic_smokenow <- cooks.distance(model.bic_smokenow); which(D_bic_smokenow > qf(0.5, 13, 400-12-1))  # 50th percentile of F with p+1 and n-p-1 df

# DFFITS
dfits_bic_smokenow <- dffits(model.bic_smokenow); which(abs(dfits_bic_smokenow) > 2*sqrt(13/400))  # 2*sqrt((p+1)/n)

# DFBETAS
dfb_bic_smokenow <- dfbetas(model.bic_smokenow); which(abs(dfb_bic_smokenow[,1]) > 2/sqrt(400))  # 2*sqrt(1/n)


# Normality check
r_bic_smokenow <- rstudent(model.bic_smokenow)
Y_hat_bic_smokenow <- predict(model.bic_smokenow)

pdf("diag_fm", height = 8, width = 16)
par(mfrow = c(2,2), mai = c(1, 0.7, 0.3, 0.3))
qqnorm(r_bic_smokenow); qqline(r_bic_smokenow)

# Standardized Residual plot for fitted y-values on bic_smokenow selected varibles, Y_hat_bic_smokenow 
plot(Y_hat_bic_smokenow, r_bic_smokenow, type = "p", xlab = "Fitted Values of BPSysAve for BIC+SmokeNow Model", ylab = "Standardized Residuals", main = "Standardized Residuals Plot with Fitted Values of BPSysAve for BIC+SmokeNow Model", col = "blue")
abline(h = 2, lty = 2); abline(h = -2, lty = 2); abline(h = 0, lwd = 0.8); lines(lowess(Y_hat_bic_smokenow, r_bic_smokenow), col="red")

# Response vs fitted values Y_hat_bic_smokenow
plot(train$BPSysAve ~ Y_hat_bic_smokenow, type="p", xlab="Fitted Values of BPSysAve for BIC+SmokeNow Model", ylab="Y, observed values of BPSysAve", main = "Scatter Plot with Y vs Fitted Values of BPSysAve for BIC+SmokeNow Model", cex.lab=1.2, col='red')
abline(lm(train$BPSysAve ~ Y_hat_bic_smokenow), lwd=2, col="blue"); lines(lowess(Y_hat_bic_smokenow, train$BPSysAve), col="green")

### Cross Validation and prediction performance of BIC + SmokeNow Model ###
ols.bic_smokenow <- ols(BPSysAve ~ ., data = train[,which(colnames(train) %in% c(sel.var.bic,"SmokeNow", "BPSysAve"))], 
                 x=T, y=T, model = T)

## 10 fold cross validation ##   
set.seed(1004079631)
bic_smokenow.cross <- calibrate(ols.bic_smokenow, method = "crossvalidation", B = 10)
## Calibration plot ##
plot(bic_smokenow.cross, las = 1, xlab = "Predicted Probability", main = "Cross-Validation calibration with BIC+SmokeNow Model")
dev.off()

## Test Error ##
pred.bic_smokenow <- predict(model.bic_smokenow, newdata = test[,which(colnames(train) %in% c(sel.var.aic,"SmokeNow", "BPSysAve"))])
## Prediction error ##
pred.error.bic_smokenow <- mean((test$BPSysAve - pred.bic_smokenow)^2)
pred.error.bic_smokenow
crit_bic_smokenow[6] = pred.error.bic_smokenow
crit_bic; crit_bic_smokenow
