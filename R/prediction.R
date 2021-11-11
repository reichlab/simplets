## functions for simple ts prediction

#' Generate predictions from a model fit with class simple_ts
#'
#' This function linearly interpolates missing values on the interior of the
#' time series, which is necessary for some model fits with moving average
#' components.  It also handles transformations and possible seasonal
#' differencing that may have been done in the fit_simple_ts function.
#'
#' @param object a model fit of class "simple_ts", as returned by fit_simple_ts
#' @param newdata numeric vector of new data to simulate forward from
#' @param horizon number of time steps forwards to simulate
#' @param nsim number of sample trajectories to simulate
#' @param seed either `NULL` or an integer that will be used in a call to
#'   `set.seed` before simulating the response vectors.  If set, the value is
#'   saved as the "seed" attribute of the returned value.  The default, `NULL`,
#'   will not change the random generator state, and return `.Random.seed`
#'   as the "seed" attribute
#' @param ... other arguments passed on to model-specific predict methods
#'
#' @return an nsim by horizon matrix with simulated values
#'
#' @export
predict.simple_ts <- function(
  simple_ts_fit,
  newdata = NULL,
  horizon = 1,
  nsim = 1,
  seed = NULL,
  ...
) {
  if(is.null(seed)) {
    seed <- .Random.seed
  } else {
    set.seed(seed)
  }
  
  transformation <- attr(simple_ts_fit, 'simple_ts_transformation')
  transform_offset <- attr(simple_ts_fit, 'simple_ts_transform_offset')
  d <- attr(simple_ts_fit, 'simple_ts_d')
  D <- attr(simple_ts_fit, 'simple_ts_D')
  ts_frequency <- attr(simple_ts_fit, 'simple_ts_ts_frequency')
  
  # Initial transformation, if necessary
  if (transformation == "box-cox") {
    bc_lambda <- attr(simple_ts_fit, 'simple_ts_bc_lambda')
  } else {
    bc_lambda <- NULL
  }

  # y from which to project forward
  if (is.null(newdata)) {
    y <- attr(simple_ts_fit, "simple_ts_y")
  } else {
    y <- newdata
  }

  transformed_y <- do_initial_transform(
    y = y,
    transformation = transformation,
    transform_offset = transform_offset,
    bc_lambda = bc_lambda)

  # Initial differencing, if necessary
  differenced_y <- do_difference(transformed_y, d = d, D = D,
    frequency = ts_frequency)
  
  # Drop leading missing values, fill in internal missing values via linear
  # interpolation.  This is necessary to ensure non-missing predictions if the
  # sarima model has a moving average component.
  # deal with internal missing or infinite values in y that can
  # result in simulated trajectories of all NAs if the model has a moving
  # average component.  Here we do this by linear interpolation.
  interpolated_y <- interpolate_and_clean_missing(differenced_y)
  
  # Sample trajectories on the transformed and differenced scale
  raw_trajectory_samples <- NextMethod(newdata = interpolated_y)
  
  # Sampled trajectories are of seasonally differenced transformed time series
  # Get to trajectories for originally observed time series ("orig") by
  # adding inverting the differencing and transformation operations
  orig_trajectory_samples <- raw_trajectory_samples
  for (i in seq_len(nsim)) {
    orig_trajectory_samples[i, ] <-
      invert_difference(
        dy = raw_trajectory_samples[i, ],
        y = transformed_y,
        d = d,
        D = D,
        frequency = ts_frequency)
    
    orig_trajectory_samples[i, ] <-
      invert_initial_transform(
        y = orig_trajectory_samples[i, ],
        transformation = transformation,
        transform_offset = transform_offset,
        bc_lambda = bc_lambda)
  }

  attr(orig_trajectory_samples, "seed") <- seed
  
  return(orig_trajectory_samples)
}
