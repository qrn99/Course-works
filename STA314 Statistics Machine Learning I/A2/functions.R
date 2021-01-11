# check existence of the required packages
needed_packages = c("tidyverse", "tidymodels")
if(!require(tidyverse)) {
  stop("The tidyverse packages must be installed. Run install.packages(\"tidyverse\") and then try again.")
}
if(!require(tidymodels)) {
  stop("The tidymodels packages must be installed. Run install.packages(\"tidymodels\") and then try again.")
}

library(tidyverse)
library(tidymodels)

#### function1 ####
fit_logistic_lasso <- function(x, y, lambda = 0, beta0 = NULL, eps = 0.0001, iter_max = 100)
{
  ## fit_logistic_lasso(x, y, lambda, beta0, eps, iter_max) takes
  ##
  ## Input:
  ## - x: matrix of predictors (not including the intercept)
  ## - y: vector of data
  ## - lambda: penalty
  ## - beta0: initial guess
  ## - eps: parameter for stopping critereon
  ## - iter_max: maximium number of iterations
  ##
  ## Output:
  ## - A list containing the members intercept, beta, and lambda
  ##
  ## Example:
  ## # make the data
  ## n = 1000
  ## dat <- tibble(x = seq(-3,3, length.out = n),
  ##               w = 3*cos(3*seq(-pi,pi, length.out = n)),
  ##               y = rbinom(n,size = 1, prob = 1/(1 + exp(-w+2*x)) )%>% as.numeric %>% factor,
  ##               cat = sample(c("a","b","c"), n, replace = TRUE)
  ## )
  ## split <- initial_split(dat, strata = c("cat"))
  ## train <- training(split)
  ## test <- testing(split)
  ## rec <- recipe(y ~ . , data = train) %>%
  ##   step_dummy(all_nominal(), -y) %>% step_zv(all_outcomes()) %>%
  ##   step_normalize(all_numeric(), -y) %>% # don't normalize y!
  ##   step_intercept() ## This is always last!
  ##
  ## ddat <- rec %>% prep(train) %>% juice
  ## x = as.matrix(select(ddat, -y, -intercept))
  ## y = ddat$y
  ##
  ## ddat_test <- rec %>% prep(test) %>% juice
  ##
  ## fit <- fit_logistic_lasso(x,y)
  ## fit
  ## $beta
  ## (intercept)           x           w       cat_b       cat_c
  ## 0.24233062 -3.73361376  2.17025090  0.09469664 -0.07367200
  ##
  ## $fct_levels
  ## [1] "0" "1"
  ##
  ## $iter
  ## [1] 7

  # try
  n <- dim(x)[1]
  # lambda = lambda*n
  x <- cbind(rep(0,n), x)  # guess intercept to be zero first
  p1 <- dim(x)[2] # first entry will be intercept, rest is beta, need this length for loop: 2:p1

  # beta0 is (0, .., 0) if not given a guess
  if (is.null(beta0)) {
    beta0 <- rep(0,p1)
  }

  ## Process the factor to be 0/1
  ## Make sure you save the names of the factor levels so we can
  ## use them in the predictions
  fct_levels <- levels(y)
  y <- as.numeric(y) - 1

  beta <- beta0

  x_beta0 <- (x %*% beta0) %>% as.numeric
  p <- exp(x_beta0)/(1 + exp(x_beta0))  ## ? change to something that does not have floating point storage issue

  # convergence loop of IRLS
  for(iter in 1:iter_max) {
    w <- p * (1 - p)
    z <- x_beta0 + (y - p)/w

    ## Compute beta with coordinate descent on penalized iteratively reweighted
    ## least squares since LASSO has 1-norm is not differentiate
    for (iter2 in 1:iter_max){ # convergence loop of coordinate descent
      for(j in 2:p1){ # inner loop of cd

        r_j = z - x[,-j] %*% beta[-j]

        # for beta[j]
        a = sum(w*(x[,j]^2))
        b = sum(w*r_j*x[,j])

        # for intercept
        c = sum(w*(z - x[,-1] %*% beta[-1]))
        d = sum(w)

        beta[j] <- sign(b)*max(0, 2*abs(b)-lambda*n)/(2*a)
        beta[1] <- c/d
      }
      # name beta
      colnames(x)[1] = c('(intercept)')
      names(beta) <- colnames(x)

      converge_cd = FALSE
      if ( max(abs(beta - beta0))/max(abs(beta)) < eps ) { #converged
        beta = beta
        # intercept = intercept
        converge_cd = TRUE
      }
      # get here means this iteration of j did not converge yet, continues
      beta0 = beta # keep track of the old beta for convergence check
      # intercept0 = intercept # keep the old intercept
    }

    # finish for looping for coordinate descent
    if (!converge_cd){
      warning(paste("Coordinate Descent Method did not converge in", iter2, "iterations", sep = " "))
    }

    ## when to stop the IRLS converge? (see below!)
    ## Compute gradient (which also computes p for the next step of the algorithm!)
    x_beta0 <- (x %*% beta) %>% as.numeric
    p <- 1/(1 + exp(-x_beta0)) # update p as intercept, beta are updated

    # compute a criterion to check whether the penalized IRLS converges or not
    # as 1-norm is not differentiate when beta = 0, we compute for different cases of beta
    check = rep(0, p1)
    for (l in 2:p1){
      if (beta[l] == 0){
        check[l] = 0
      } else {
        check[l] = sign(beta[l])*n*lambda  # -lambda if beta < 0, +lambda if beta > 0
      }
    }

    grad <- -2*t(x) %*% (y - p) + check

    if (sqrt(sum(grad^2)) < eps) {
        return(list(beta = beta, fct_levels = fct_levels, iter = iter))
      }
      ## Otherwise we go again!
      ## We don't need to do anything here - we've already computed x_beta0 and p
      ## for the next iteration!
  }

  warning(paste("IRLS Method did not converge in", iter_max, "iterations", sep = " "))
    return(list(beta = beta, fct_levels = fct_levels, iter = iter))
}

#### function2 ####
######################Is this how you write predict...
predict_logistic_lasso <- function(fit, new_x)
{
  ## predict_logistic_lasso(fit, new_x) takes
  ##
  ## Input:
  ## - fit: Output from fit_logistic_lasso
  ## - new_x: Data to predict at (may be more than one point!)
  ##
  ## Output:
  ## - A numeric vector that is the data being predicted based on the input fit for beta and intercept
  ##
  ## Example:
  ## # continue from previous function example
  ## > fit <- fit_logistic_lasso(x,y)
  ## > predict <- predict_logistic_lasso(fit, new_x)
  ## [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
  ## [48] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 1 1 1
  ## [95] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
  ## [142] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  ## [189] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  ## [236] 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  ## Levels: 0 1
  ##
  intercept <- fit$beta[1]
  beta_coef <- fit$beta[-1]

  numeric_pred <- (new_x %*% beta_coef >= 0) %>% as.numeric
  return(fit$fct_levels[numeric_pred + 1] %>% factor)

}