#' Compute cumulative excess deaths
#' 
#' This function takes the output of the `excess_model` function, a start date, and 
#' end date and calculates excess mortality and standard errors.
#' 
#' @param fit The output of `excess_model` 
#' @param start The start date 
#' @param end The end date
#' 
#' @return A data frame with the following columns
#' \describe{
#' \item{date}{The date}
#' \item{observed}{The observed deficit count,which is the sum of observed minus expected until that date}
#' \item{sd}{The standard deviation for the deficit count for that date if year is typical}
#' \item{fitted}{The estimated of edeficit count based on the estimated smooth event effect curve}
#' \item{se}{The standard error for `fitted`}
#' }
#'
#' @examples
#' 
#' data(new_jersey_counts)
#' exclude_dates <- as.Date("2012-10-29") + 0:180
#' control_dates <- seq(min(new_jersey_counts$date), min(exclude_dates) - 1, by = "day")
#' f <- excess_model(new_jersey_counts, 
#' start = as.Date("2012-09-01"), 
#' end = as.Date("2013-09-01"),
#' exclude = exclude_dates,
#' model = "correlated",
#' weekday.effect = TRUE,
#' control.dates = control_dates)
#' 
#' excess_cumulative(f, 
#' start = as.Date("2017-12-15"), 
#' end = as.Date("2017-12-21") )
#' 
#' @export
deficit_cumulative <- function(fit, start, end){
  if (!"curve_fit" %in% attr(fit, "type"))
    stop("This is not the correct deficit_model fit, needs curve fit.")

  ind             <- which(fit$date %in% seq(start, end, by = "day"))
  n               <- length(ind) # returns the # of elements in the ind object
  A               <- matrix(1, n, n)
  A[upper.tri(A)] <- 0                                         # It sets all the upper triangular elements of matrix A (excluding the diagonal) to 0
  A               <- sweep(A, 2, fit$expected[ind], FUN = "*") # This line multiplies each column (2 means columns) of the matrix A by a corresponding value from fit$expected[ind]
                                                               # So, now it's scaling each column of A by expected values from the model fit (fit$expected[ind])
  fhat <- matrix(fit$fitted[ind], ncol = 1)                    # creates a column vector (matrix with 1 column) called fhat from the values in fit$fitted[ind]

 
   # Below, this code is constructing a data frame with:
  
  fit_deficit   <- A %*% fhat                                    # multiply design matrix A by fitted values fhat to get cumulative modeled effect.
  obs_deficit   <- cumsum(fit$observed[ind] - fit$expected[ind])
  varcov       <- fit$x[ind,] %*% fit$betacov %*% t(fit$x[ind,]) # Compute variance-covariance matrix of the linear predictor based on model design matrix x and estimated betacov
  diag(varcov) <- fit$se[ind]^2                                  # Replace diagonal entries of varcov with squared standard errors — probably to correct for over/under-dispersion or model misspecification.
  fitted_se    <- sqrt(diag(A %*%  varcov %*% t(A)))
  sd           <- sqrt(diag(A %*% (fit$cov[ind,ind] + diag(fit$log_expected_se[ind]^2)) %*% t(A))) # this computes the variance of the observed cumulative deficit.
                                                                                                   # It adds fit$cov (covariance of observed data) and extra uncertainty from fit$log_expected_se.
  data.frame(date     = fit$date[ind], # The time points
             observed = obs_deficit,   # Cumulative observed - expected difference
             sd       = sd,            # Modeled cumulative deficit/excess from the model.
             fitted   = fit_deficit,   # Standard deviation of the observed cumulative deficit.
             se       = fitted_se)     # Standard error of the fitted cumulative deficit.
}
