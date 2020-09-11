rm(list = ls())
setwd("~/Desktop/UofT/19_20/STA302")

## Q1 ##
## ---- setup1
set.seed(1004079631)
beta0 <- rnorm(n = 1, mean = 0, sd = 1)  # The population beta_0
beta1 <- runif(n = 1, min = 1, max = 3)  # The population beta_1
sig2 <- rchisq(n = 1, df = 25)  # The error variance sigma^2

nsample <- 5  # Sample size
n.sim <- 1000  # The number of simulations
sigX <- 0.2  # The variances of X

## ---- setup2
X <- rnorm(n = nsample, mean = 0, sd = sqrt(sigX))  #Simulate the predictor variable

b0 <- vector()  # saves the sample estimates of beta_0
b1 <- vector()  # saves the sample estimates of beta_1
sig2hat <- vector()  # saves the sample estimates of error variance sigma^2

for(i in 1:n.sim){  # Assign the estimators for $n.sim$ samples individually
  Y <- beta0 + beta1*X + rnorm(n = nsample, mean = 0, sd = sqrt(sig2))
  model <- lm(Y ~ X)
  b0[i] <- coef(model)[1]
  b1[i] <- coef(model)[2]
  sig2hat[i] <- summary(model)$sigma^2
}


## ---- mean_of_esi_1000
## Q2 ##
# The mean of esitmators b0, b1, and error variance s^2 for sample size = nsample with 100 simulations
mean(b0); mean(b1); mean(sig2hat)
# The true population parameters beta_0, beta_1, and error variance sigma^2
beta0; beta1; sig2


## ---- hiso_b0b1s2
## Q3 ##
# Construct hisotgrams of estiamtors for b0, b1, and error variance s^2 with $n.sim$ samples
par(mar=c(1,1,1,1)); par(mfrow = c(3,3))  # set up the margin of the plots
hist(b0); hist(b1); hist(sig2hat)

## ---- var_of_reg_paras
## Q4 ##
# Calculate the true variance of b0 and b1
sumx <- sum(X); xbar <- sumx/nsample
sumx2 <- sum(X^2); SXX <- sumx2 - nsample*(xbar^2)
var_beta_0 <- sig2*(1/nsample + xbar^2/SXX)
var_beta_1 <- sig2/SXX

var_beta_0; var_beta_1  # true variance of b0 and b1

# Save sample variance of b0, b1 for each simulation
var_b0 <- vector()
var_b1 <- vector()
for(i in 1:n.sim){  # Assign the var estimators of b0, b1 for $n.sim$ samples individually
  var_b0[i] <- sig2hat[i]*(1/nsample + xbar^2/SXX)
  var_b1[i] <- sig2hat[i]/SXX
}

mean(var_b0); mean(var_b1)  # the mean of the sample variances of b0 and b1

## ---- 95%_CI_z
## Q5 ##
# 95% Z-test CI, call it CI_z
ll_b0_z <- vector(); ul_b0_z <- vector()  # Save lower level, upper level of CI_z for b0
ll_b1_z <- vector(); ul_b1_z <- vector()  # Save lower level, upper level of CI_z for b1
counts_b0_z <- 0; counts_b1_z <- 0  # Save how many CI_z the true reg parameters do fall into
z_95 <- qnorm(.025, lower.tail=FALSE)  # z score for 95% CI

for(i in 1:n.sim){  # calculate 95% CI (Z-test) for b0, b1
  ll_b0_z[i] <- b0[i] - z_95*sqrt(var_beta_0); ul_b0_z[i] <- b0[i] + z_95*sqrt(var_beta_0)
  ll_b1_z[i] <- b1[i] - z_95*sqrt(var_beta_1); ul_b1_z[i] <- b1[i] + z_95*sqrt(var_beta_1)

  # Check if CI_z contains the true value of the parameters
  if ((ll_b0_z[i] <= beta0) && (beta0 <= ul_b0_z[i])) {
    counts_b0_z <- counts_b0_z + 1
  }
  if ((ll_b1_z[i] <= beta1) && (beta1 <= ul_b1_z[i])) {
    counts_b1_z <- counts_b1_z + 1
  }
}

counts_b0_z/n.sim; counts_b1_z/n.sim # Precentatge of the true value of regression parameters fall into CI_z

## ---- 95%_CI_t
# 95% T-test CI, call it CI_t
# Save lower level, uppler level of CI_t for b0, b1
ll_b0_t <- vector(); ul_b0_t <- vector()
ll_b1_t <- vector(); ul_b1_t <- vector()
counts_b0_t <- 0; counts_b1_t <- 0  # Save how many CI_t the true reg parameters do fall into
t_95 <- qt(c(.025, .975), df=nsample-2)[2] # t score for 95% CI with sample size = nsample = 5, so df = nsample-2

for(i in 1:n.sim){  #calculate 95% CI (T-test) for b0, b1
  # b0[i] +- t_95 * var_b0[i] is the CI_t for b0[i], as var_b0[i] is the sample var calculated before for b0
  ll_b0_t[i] <- b0[i] - t_95 * sqrt(var_b0[i]); ul_b0_t[i] <- b0[i] + t_95 * sqrt(var_b0[i])
  # b1[i] +- t_95 * var_b1[i] is the CI_t for b1[i], as var_b1[i] is the sample var calculated before for b1
  ll_b1_t[i] <- b1[i] - t_95 * sqrt(var_b1[i]); ul_b1_t[i] <- b1[i] + t_95 * sqrt(var_b1[i])
  
  # Check if CI_t contains the true value of the parameters
  if ((ll_b0_t[i] <= beta0) && (beta0 <= ul_b0_t[i])) {
    counts_b0_t <- counts_b0_t + 1
  }
  if ((ll_b1_t[i] <= beta1) && (beta1 <= ul_b1_t[i])) {
    counts_b1_t <- counts_b1_t + 1
  }
}

counts_b0_t/n.sim; counts_b1_t/n.sim # Precentatge of the true value of regression parameters fall into CI_t

## ---- Increase_samp_size_25
## Q6 ##
# Start with sample size 25
set.seed(1004079631)
X_25 <- rnorm(n = 25, mean = 0, sd = sqrt(sigX))  #Simulate the predictor variable with new sample size
b0_25 <- vector()  # saves the sample estimates of beta_0 with sample size 25
b1_25 <- vector()  # saves the sample estimates of beta_1 with sample size 25
sig2hat_25 <- vector()  # saves the sample estimates of error variance sigma^2 with sample size 25

for(i in 1:n.sim){  # Assign the estimators for $n.sim$ samples individually
  Y_25 <- beta0 + beta1*X_25 + rnorm(n = 25, mean = 0, sd = sqrt(sig2))
  model_25 <- lm(Y_25 ~ X_25)
  b0_25[i] <- coef(model_25)[1]
  b1_25[i] <- coef(model_25)[2]
  sig2hat_25[i] <- summary(model_25)$sigma^2
}

## ---- Increase_samp_size_50
# Increase sample size to 50
set.seed(1004079631)
X_50 <- rnorm(n = 50, mean = 0, sd = sqrt(sigX))  #Simulate the predictor variable with new sample size
b0_50 <- vector()  # saves the sample estimates of beta_0 with sample size 50
b1_50 <- vector()  # saves the sample estimates of beta_1 with sample size 50
sig2hat_50 <- vector()  # saves the sample estimates of error variance sigma^2 with sample size 50

for(i in 1:n.sim){  # Assign the estimators for $n.sim$ samples individually with sample size 50
  Y_50 <- beta0 + beta1*X_50 + rnorm(n = 50, mean = 0, sd = sqrt(sig2))
  model_50 <- lm(Y_50 ~ X_50)
  b0_50[i] <- coef(model_50)[1]
  b1_50[i] <- coef(model_50)[2]
  sig2hat_50[i] <- summary(model_50)$sigma^2
}

## ---- Increase_samp_size_100
# Increase sample size to 100
set.seed(1004079631)
X_100 <- rnorm(n = 100, mean = 0, sd = sqrt(sigX))  #Simulate the predictor variable with new sample size
b0_100 <- vector()  # saves the sample estimates of beta_0 with sample size 100
b1_100 <- vector()  # saves the sample estimates of beta_1 with sample size 100
sig2hat_100 <- vector()  # saves the sample estimates of error variance sigma^2 with sample size 100

for(i in 1:n.sim){  # Assign the estimators for $n.sim$ samples individually with sample size 100
  Y_100 <- beta0 + beta1*X_100 + rnorm(n = 100, mean = 0, sd = sqrt(sig2))
  model_100 <- lm(Y_100 ~ X_100)
  b0_100[i] <- coef(model_100)[1]
  b1_100[i] <- coef(model_100)[2]
  sig2hat_100[i] <- summary(model_100)$sigma^2
}

## ---- mean_of_esi_1000_samp_size_inc
# The mean of esitmators b0, b1, and error variance s^2 for sample size = 25 with 100 simulations
mean(b0_25); mean(b1_25); mean(sig2hat_25)
# The mean of esitmators b0, b1, and error variance s^2 for sample size = 50 with 100 simulations
mean(b0_50); mean(b1_50); mean(sig2hat_50)
# The mean of esitmators b0, b1, and error variance s^2 for sample size = 100 with 100 simulations
mean(b0_100); mean(b1_100); mean(sig2hat_100)
# The true population parameters beta_0, beta_1, and error variance sigma^2
beta0; beta1; sig2

## ---- hiso_b0b1s2_samp_size_inc
# Construct hisotgrams of sample estiamtors to b0, b1, and error variance sig2 with increasing sample size
par(mar=c(1,1,1,1)); par(mfrow = c(4,3), mai = c(0.3, 0.1, 0.1, 0.1))  # Side by side comparsion
hist(b0_25); hist(b1_25); hist(sig2hat_25)
hist(b0_50); hist(b1_50); hist(sig2hat_50)
hist(b0_100); hist(b1_100); hist(sig2hat_100)

## ---- var_of_reg_paras_samp_size_inc_25
# Calculate the true variance of b0 and b1 for sample size 25
sumx_25 <- sum(X_25); xbar_25 <- sumx_25/25
sumx2_25 <- sum(X_25^2); SXX_25 <- sumx2_25 - 25*xbar_25^2
var_beta_0_25 <- sig2*(1/25 + xbar_25^2/SXX_25)
var_beta_1_25 <- sig2/SXX_25

var_beta_0_25; var_beta_1_25  # true variance of b0 and b1 for sample size 25

# Calculate sample variance of b0, b1 for each simulation for sample size 25
var_b0_25 <- vector()
var_b1_25 <- vector()
for(i in 1:n.sim){  # Assign the var estimators of b0, b1 for $n.sim$ samples individually
  var_b0_25[i] <- sig2hat_25[i]*(1/25 + xbar_25^2/SXX_25)
  var_b1_25[i] <- sig2hat_25[i]/SXX_25
}

mean(var_b0_25); mean(var_b1_25)  # the mean of sample variances of  b0 and b1

## ---- var_of_reg_paras_samp_size_inc_50
# Calculate the true variance of b0 and b1 for sample size 50
sumx_50 <- sum(X_50); xbar_50 <- sumx_50/50
sumx2_50 <- sum(X_50^2); SXX_50 <- sumx2_50 - 50*xbar_50^2
var_beta_0_50 <- sig2*(1/50 + xbar_50^2/SXX_50)
var_beta_1_50 <- sig2/SXX_50

var_beta_0_50; var_beta_1_50  # true variance of b0 and b1 for sample size 50

# Calculate sample variance of b0, b1 for each simulation for sample size 50
var_b0_50 <- vector()
var_b1_50 <- vector()
for(i in 1:n.sim){  # Assign the var estimators of b0, b1 for $n.sim$ samples individually
  var_b0_50[i] <- sig2hat_50[i]*(1/50 + xbar_50^2/SXX_50)
  var_b1_50[i] <- sig2hat_50[i]/SXX_50
}

mean(var_b0_50); mean(var_b1_50)  # the mean of sample variances of b0 and b1

## ---- var_of_reg_paras_samp_size_inc_100
# Calculate the true variance of b0 and b1 for sample size 100
sumx_100 <- sum(X_100); xbar_100 <- sumx_100/100
sumx2_100 <- sum(X_100^2); SXX_100 <- sumx2_100 - 100*xbar_100^2
var_beta_0_100 <- sig2*(1/100 + xbar_100^2/SXX_100)
var_beta_1_100 <- sig2/SXX_100

var_beta_0_100; var_beta_1_100  # true variance of b0 and b1 for sample size 100

# Calculate sample variance of b0, b1 for each simulation for sample size 100
var_b0_100 <- vector()
var_b1_100 <- vector()
for(i in 1:n.sim){  # Assign the var estimators of b0, b1 for $n.sim$ samples individually
  var_b0_100[i] <- sig2hat_100[i]*(1/100 + xbar_100^2/SXX_100)
  var_b1_100[i] <- sig2hat_100[i]/SXX_100
}

mean(var_b0_100); mean(var_b1_100)  # the mean of sample variances of b0 and b1 for sample size 100


## ---- change_error_var_sig2
## Q7 ##
set.seed(1004079631)
sig2_small <- 0.012
X_small <- rnorm(n = 100, mean = 0, sd = sqrt(sigX))  #Simulate the predictor variable with sample size = 100

b0_small <- vector()  # saves the sample estimates of beta_0 with small error var
b1_small <- vector()  # saves the sample estimates of beta_1 with small error var
sig2hat_small <- vector()  # saves the sample estimates of new error variance sig2_small


for(i in 1:n.sim){  # Assign the estimators for $n.sim$ samples individually with new error variance sig2_small
  Y_small <- beta0 + beta1*X_small + rnorm(n = 100, mean = 0, sd = sqrt(sig2_small))
  model_small <- lm(Y_small ~ X_small)
  b0_small[i] <- coef(model_small)[1]
  b1_small[i] <- coef(model_small)[2]
  sig2hat_small[i] <- summary(model_small)$sigma^2
}

# Increase error var to large
set.seed(1004079631)
sig2_large <- 1100.289
X_large <- rnorm(n = 100, mean = 0, sd = sqrt(sigX))  #Simulate the predictor variable with sample size = 100

b0_large <- vector()  # saves the sample estimates of beta_0 with large error var
b1_large <- vector()  # saves the sample estimates of beta_1 with large error var
sig2hat_large <- vector()  # saves the sample estimates of new error variance sig2_large

for(i in 1:n.sim){  # Assign the estimators with new error variance sig2_large
  Y_large <- beta0 + beta1*X_large + rnorm(n = 100, mean = 0, sd = sqrt(sig2_large))
  model_large <- lm(Y_large ~ X_large)
  b0_large[i] <- coef(model_large)[1]
  b1_large[i] <- coef(model_large)[2]
  sig2hat_large[i] <- summary(model_large)$sigma^2
}

## ---- mean_of_esi_1000_error_var_inc
# The mean of esitmators b0, b1 for small error var
mean(b0_small); mean(b1_small)
# The mean of esitmators b0, b1 for large error var
mean(b0_large); mean(b1_large)
# The true population parameters beta_0, beta_1
beta0; beta1

## ---- hiso_b0b1s2_error_var_inc
# Construct hisotgrams of estiamtors for b0 & b1, small error var vs large error var
par(mar=c(1,1,1,1)); par(mfrow = c(3,2), mai = c(0.3, 0.1, 0.1, 0.1))  # set up the marigins
hist(b0_small); hist(b1_small)  # Side by Side comparsion
hist(b0_large); hist(b1_large)

## ---- var_of_reg_paras_error_var_inc_small
# Calculate the true variance of b0 and b1 for small error var
sumx_small <- sum(X_small); xbar_small <- sumx_small/100
sumx2_small <- sum(X_small^2); SXX_small <- sumx2_small - 100*xbar_small^2
var_beta_0_small <- sig2*(1/100 + xbar_small^2/SXX_small)
var_beta_1_small <- sig2/SXX_small

var_beta_0_small; var_beta_1_small  # true variance of b0 and b1 for small error var

# Calculate sample variance of b0, b1 for for small error var
var_b0_small <- vector()
var_b1_small <- vector()
for(i in 1:n.sim){  # Assign the var estimators of b0, b1 for $n.sim$ samples individually
  var_b0_small[i] <- sig2hat_small[i]*(1/100 + xbar_small^2/SXX_small)
  var_b1_small[i] <- sig2hat_small[i]/SXX_small
}

mean(var_b0_small); mean(var_b1_small)  # the mean of sample variances of b0 and b1 for small error var

## ---- var_of_reg_paras_error_var_inc_large
# Calculate the true variance of b0 and b1 for large error var
sumx_large <- sum(X_large); xbar_large <- sumx_large/100
sumx2_large <- sum(X_large^2); SXX_large <- sumx2_large - 100*xbar_large^2
var_beta_0_large <- sig2*(1/100 + xbar_large^2/SXX_large)
var_beta_1_large <- sig2/SXX_large

var_beta_0_large; var_beta_1_large  # true variance of b0 and b1 for large error var

# Calculate sample variance of b0, b1 for for large error var
var_b0_large <- vector()
var_b1_large <- vector()
for(i in 1:n.sim){  # Assign the var estimators of b0, b1 for $n.sim$ samples individually
  var_b0_large[i] <- sig2hat_large[i]*(1/100 + xbar_large^2/SXX_large)
  var_b1_large[i] <- sig2hat_large[i]/SXX_large
}

mean(var_b0_large); mean(var_b1_large)  # the mean of sample variances of b0 and b1 for large error var

