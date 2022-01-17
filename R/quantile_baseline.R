#' Check if object is of class quantile_baseline
#'
#' @param object an object that may be a quantile_baseline object
#'
#' @return boolean; whether object is inherits quantile_baseline class
#'
#' @export
is.quantile_baseline <- function(object) {
  return(inherits(object, "quantile_baseline"))
}


#' Create a quantile_baseline object
#'
#' @param inc_diffs historical first differences in incidence
#' @param symmetrize logical. if TRUE (the default), we collect the first
#' differences of incidence and their negatives; the resulting distribution on
#' differences is symmetric. If FALSE, we use only the observed differences.
#'
#' @return quantile_baseline fit object
#'
#' @export
new_quantile_baseline <- function(inc_diffs, symmetrize = TRUE) {
  if(symmetrize) {
    quantile_baseline <- structure(
      c(inc_diffs, -inc_diffs),
      symmetrize = symmetrize,
      class = 'quantile_baseline'
    )
  } else {
    quantile_baseline <- structure(
      inc_diffs,
      symmetrize = symmetrize,
      class = 'quantile_baseline'
    )
  }
  
  return(quantile_baseline)
}


#' Fit a quantile baseline model to historical disease incidence
#'
#' @param incidence numeric vector of disease incidence in past time points
#' @param symmetrize logical. if TRUE (the default), we collect the first
#' differences of incidence and their negatives; the resulting distribution on
#' differences is symmetric. If FALSE, we use only the observed differences.
#' @param window_size integer optional number of past time points to use for
#'   finding first differences.  If not provided, all past first differences
#'   will be used.
#' @param ... other arguments are ignored
#'
#' @return quantile_baseline fit object
#'
#' @export
fit_quantile_baseline <- function(
    incidence,
    symmetrize = TRUE,
    window_size = length(incidence) - 1,
    ...) {
  if(is.na(window_size)) {
    window_size <- length(incidence) - 1
  }
  if(window_size >= length(incidence)) {
    window_size <- length(incidence) - 1
  }
  diffs <- tail(diff(incidence), window_size)
  diffs <- diffs[!is.na(diffs)]
  return(new_quantile_baseline(
    inc_diffs=diffs,
    symmetrize=symmetrize))
}


#' Predict future disease incidence by resampling one-step-ahead forecasts
#'
#' @param quantile_baseline a quantile_baseline fit object
#' @param newdata numeric vector of length at least one with incident counts
#' @param cum_data numeric vector of length at least one with cumulative counts
#' @param quantiles quantile levels for which  to generate predictions
#' @param horizon number of time steps forward to predict
#' @param nsim number of samples to use for generating predictions at
#' horizons greater than 1
#' @param origin string specifying whether to project forward from the
#' most recent observation (`"obs"`) or from the fitted value from a LOESS
#' smooth (`"loess"`)
#'
#' @return matrix of samples of incidence
#'
#' @importFrom tsibble as_tsibble
#' @importFrom fabletools model
#' @importFrom feasts STL
#'
#' @export
predict.quantile_baseline <- function(
  quantile_baseline,
  newdata,
  quantiles,
  horizon,
  nsim,
  origin = c("obs", "loess")
  ) {
  origin <- match.arg(origin)
  
  symmetrize <- attr(quantile_baseline, "symmetrize")

  results <- matrix(NA, nrow = nsim, ncol = horizon)
  
  if (origin == "obs") {
    pred_origin <- tail(newdata, 1)
  } else {
    stl_formula <- y ~ trend(window = 7) +
      season(period = 1, window = 1)

    pred_origin <- data.frame(
        y = newdata,
        time = Sys.Date() - rev(seq_along(newdata))) %>%
      as_tsibble(index = time) %>%
      model(STL(stl_formula, robust = TRUE)) %>%
      generics::components() %>%
      as_tibble() %>%
      pull(trend) %>%
      tail(1)
  }
  
  # Case for horizon 1 is different because sampling is not necessary; we can
  # extract exact quantiles
  
  ## sample incidence, then correct it:
  ## - enforce median difference is 0
  ## - enforce incidence is non-negative
  sampled_inc_diffs <- quantile(
    quantile_baseline,
    probs = seq(from = 0, to = 1.0, length = nsim))
  sampled_inc_raw <- pred_origin + sampled_inc_diffs
  
  ## save as a column in results
  results[, 1] <- sampled_inc_raw
  
  for (h in (1 + seq_len(horizon - 1))) {
    sampled_inc_diffs <- sample(sampled_inc_diffs, size = nsim, replace = FALSE)
    sampled_inc_raw <- sampled_inc_raw + sampled_inc_diffs
    
    # force median difference = 0
    if (symmetrize) {
      sampled_inc_corrected <- sampled_inc_raw -
        (median(sampled_inc_raw) - pred_origin)
    } else {
      sampled_inc_corrected <- sampled_inc_raw
    }
    
    ## save as a column in results
    results[, h] <- sampled_inc_corrected
  }
  
  return(results)
}
