## functions for estimating simple ts models

#' Estimate model
#'
#' @param y a univariate time series or numeric vector.
#' @param ts_frequency frequency of time series.  Must be provided if y is not
#'   of class "ts".  See the help for stats::ts for more.
#' @param model string specifying model to fit: one of 'quantile_baseline' or
#'   'spline_smoother'
#' @param transformation character specifying transformation type:
#'   "box-cox", "sqrt", "log", or "none".  See details for more.
#' @param transform_offset numeric offset used before the Box-Cox and log
#'   transformations; the offset is added to all observations before
#'   transforming.  Default value of 0.5 allows us to use the Box-Cox and log
#'   transforms (which require positive inputs) in case of observations of 0,
#'   and also ensures that the reverse transformed values will always be at
#'   least -0.5, so that they round up to non-negative values.
#' @param d integer order of first differencing; default is 0
#' @param D integer order of seasonal differencing; default is 0
#' @param ... arguments passed on to model-specific fit method
#'
#' @return a simple_ts model fit
#'
#' @details This function is a wrapper around model-specific fit methods,
#'   providing some preliminary transformations of the data.
#'   Formal and informal experimentation has shown these preliminary
#'   transformations to be helpful with a few infectious disease time series
#'   data sets.
#'
#' @importFrom car powerTransform
#'
#' @export
fit_simple_ts <- function(
  y,
  ts_frequency = 1,
  model = 'quantile_baseline',
  transformation = 'box-cox',
  transform_offset = 0.5,
  d = 0,
  D = 0,
  ...) {
  # Validate arguments
  if (!(is.numeric(y) || is.ts(y))) {
    stop("The argument y must be a numeric vector or object of class 'ts'.")
  }
  
  if (is.ts(y)) {
    ts_frequency <- frequency(y)
  }
  
  model <- match.arg(model, c('quantile_baseline'))
  transformation <- match.arg(transformation, c('none', 'log', 'sqrt', 'box-cox'))
  
  # Initial transformation, if necessary
  if (identical(transformation, "box-cox")) {
    est_bc_params <- powerTransform(
      y + transform_offset,
      family = "bcPower")
    est_bc_lambda <- est_bc_params$lambda
  }
  transformed_y <- do_initial_transform(
    y = y,
    transformation = transformation,
    transform_offset = transform_offset,
    bc_lambda = est_bc_lambda)

  # Initial differencing, if necessary
  differenced_y <- do_difference(transformed_y, d = d, D = D,
    frequency = ts_frequency)
  
  # Get fit
  if (model == 'quantile_baseline') {
    simple_ts_fit <- fit_quantile_baseline(incidence = differenced_y, ...)
  }
  
  # Save information needed for prediction
  attr(simple_ts_fit, 'simple_ts_call') <- match.call()
  params_to_save <- c("y", "ts_frequency", "transformation", "transform_offset",
                      "d", "D")
  for (param_name in params_to_save) {
    attr(simple_ts_fit, paste0("simple_ts_", param_name)) <- get(param_name)
  }
  if (transformation == "box-cox") {
    attr(simple_ts_fit, 'simple_ts_bc_lambda') <- est_bc_lambda
  }
  
  class(simple_ts_fit) <- c("simple_ts", class(simple_ts_fit))
  
  return(simple_ts_fit)
}
