# set directory to have access to functions.R on my local computer
# setwd("/Users/nancyqiu/OneDrive - University of Toronto/UofT/20_21/STA314/a2")

# use functions in functions.R
source("functions.R")

logistic_lasso <- function(mode = "classification", penalty) 
{
  ## logistic_LASSO(mode = "classification", penalty) gets logistic_lasso generated in functions.R 
  ## into a form where we can use it with tidymodels and parsnip model.
  ## 
  ## Input:
  ## - mode: the mode of the model, default "classification" since we are doing logistic regression
  ## - penalty: the lambda value the penalize the logistic lasso regression
  ##
  ## Output:
  ## – void, just prepares the model into tidymodel package as logistic_lasso() %>% set_engine("fit_logistic_lasso") 
  ##
  ## Example:
  ## > logistic_lasso(mode = "classification", penalty)
  ## 
  args <- list(penalty = rlang::enquo(penalty))
  new_model_spec("logistic_lasso",
                 args = args,
                 mode = mode,
                 eng_args = NULL,
                 method = NULL,
                 engine = NULL)
}

# Start setting up the model
set_new_model("logistic_lasso")
set_model_mode(model = "logistic_lasso", mode = "classification")

set_model_engine("logistic_lasso",
                 mode = "classification",
                 eng = "fit_logistic_lasso"
)

set_dependency("logistic_lasso", eng = "fit_logistic_lasso", pkg = "base")

set_model_arg(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  parsnip = "penalty", ## what parsnip will call it
  original = "lambda", ## what we call it!
  func = list(pkg = "dials", fun = "penalty"), ## Use dials::penalty() to set
  has_submodel = FALSE # If you don't know, don't worry.
)

set_encoding(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  mode = "classification",
  options = list(
    predictor_indicators = "traditional",
    compute_intercept = TRUE,  # We are cd for intercept as we are penalizing it
    remove_intercept = TRUE, # x do not have intercept
    allow_sparse_x = FALSE
  )
)


# set up fitting
set_fit(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  mode = "classification",
  value = list(
    interface = "matrix",
    protect = c("x", "y"),
    func = c(fun = "fit_logistic_lasso"),
    defaults = list()
  )
)

# set up prediction
set_pred(
  model = "logistic_lasso",
  eng = "fit_logistic_lasso",
  mode = "classification",
  type = "class", #different than ridge
  value = list(
    pre = NULL,
    post = NULL,
    func = c(fun = "predict_logistic_lasso"),
    args = list(
      fit = expr(object$fit),
      new_x = expr(as.matrix(new_data[, names(object$fit$beta[-1])]))
    )
  )
)

# finalize the model
update.logistic_lasso <- function(object, penalty = NULL, ...) 
{
  ## update.logistic_lasso() finalizes the model workflow (will be used)
  ## 
  ## Input:
  ## - object: the output of fit_logistic_lasso
  ##
  ## Output:
  ## – void, just finalizes the model
  ##
  ## Example:
  ## > update.logistic_lasso(fit_logistic_lasso(...))
  if(! is.null(penalty)) {
    object$args <- list(penalty = enquo(penalty))
  }
  new_model_spec("logistic_lasso", args = object$args, eng_args = NULL,
                 mode = "classification", method = NULL, engine = object$engine)
}


