#  Adapting the excess_model function for analyzing deficit hospitalizations:
 
excess_model <- function(counts,
                         start = NULL,
                         end = NULL,
                         knots.per.year = 12,
                         event = NULL,
                         intervals = NULL,
                         discontinuity = TRUE,
                         model = c("quasipoisson", "poisson", "correlated"),
                         exclude = NULL,
                         include.trend = TRUE,
                         trend.knots.per.year = 1/7,
                         harmonics = 2,
                         frequency = NULL,
                         weekday.effect = FALSE,
                         control.dates = NULL,
                         max.control = 5000,
                         order.max = 14,
                         aic = TRUE,
                         maxit = 25,
                         epsilon = 1e-8,
                         alpha = 0.05,
                         min.rate = 0.0001,
                         keep.counts = FALSE,
                         keep.components = TRUE,
                         verbose = TRUE) {
  
  if (!"compute_expected" %in% class(counts)) {
    if (verbose) message("Computing expected counts.")
    counts <- compute_expected(counts,
                               exclude = exclude,
                               include.trend = include.trend,
                               trend.knots.per.year = trend.knots.per.year,
                               harmonics = harmonics,
                               frequency = frequency,
                               weekday.effect = weekday.effect,
                               keep.components = keep.components,
                               verbose = verbose)
  }
  
  correlated.errors <- match.arg(model) == "correlated"
  frequency <- attr(counts, "frequency")
  dispersion <- attr(counts, "dispersion")
  
  if (frequency == 12 & correlated.errors) {
    stop("Correlated error model cannot be fitted with monthly data.")
  }
  
  if (any(counts$excluded)) {
    exclude <- counts$date[counts$excluded] 
  } else {
    exclude <- NULL
  }
  
  if (correlated.errors & is.null(control.dates)) {
    if (!is.null(exclude)) {
      warning("No control region supplied. Using data up to the first excluded point as control.")
      control.dates <- seq(min(counts$date), min(exclude, na.rm = TRUE) - 1, by = "day")
    } else {
      warning("No control region supplied. Using all available data as control.")
      control.dates <- counts$date
    }
  }
  
  if ((is.null(start) & !is.null(end)) | (!is.null(start) & is.null(end))) {
    stop("Both start and end dates must be provided.")
  }
  
  if (correlated.errors) {
    arfit <- fit_ar(counts, control.dates, order.max = order.max, aic = aic)
    s <- arfit$sigma
    if (verbose) message("AR model order selected: ", length(arfit$ar), ", Residual standard error: ", signif(arfit$sigma, 2))
  }
  
  if (!is.null(start) & !is.null(end)) {
    include_dates = seq(start, end, by = "day")
    if (frequency == 365)
      include_dates <- include_dates[!(lubridate::month(include_dates) == 2 & lubridate::day(include_dates) == 29)]
    
    ind <- which(counts$date %in% include_dates)
    date <- counts$date[ind]
    n <- length(ind)
    x <- (0:(n - 1)) / frequency * 365
    mu <- pmax(min.rate, counts$expected[ind])
    obs <- counts$outcome[ind]
    pop <- counts$population[ind]
    
    nknots <- round(knots.per.year * as.numeric(max(date) - min(date)) / 365)
    knots <- x[round(seq(1, n, length = nknots + 2))]
    knots <- knots[-c(1, length(knots))]
    
    X <- cbind(1, splines::ns(x, knots = knots))
    
    fit <- glm(obs ~ X - 1, offset = log(mu), family = quasipoisson(link = "log"))
    tmp <- predict(fit, se.fit = TRUE, type = "response")
    fhat <- pmax(tmp$fit / mu - 1, min.fhat)
    lambda <- fitted.values(fit)
    cova <- summary(fit)$cov.unscaled * dispersion
    lambda_vari <- lambda^2 * diag(X %*% cova %*% t(X))
    mu_vari <- counts$expected[ind]^2 * counts$log_expected_se[ind]^2
    se <- sqrt((lambda_vari / mu^2) + (lambda^2 * mu_vari / mu^4))
    Sigma <- diag(n) * summary(fit)$dispersion / mu
    
    ## Define the thresholds for detecting significant deviations
    excess_threshold <- qnorm(1 - alpha / 2) * se
    deficit_threshold <- -excess_threshold
    
    ## Identify intervals of excess and deficit
    excess_ind <- which(fhat - excess_threshold >= 0)
    deficit_ind <- which(fhat + deficit_threshold <= 0)
    
    ## Process identified intervals
    excess_intervals <- identify_intervals(excess_ind, date)
    deficit_intervals := identify_intervals(deficit_ind, date)
    
    ret <- list(
      date = date,
      observed = obs,
      expected = mu,
      fitted = mu * (1 + fhat),
      se = se,
      excess_intervals = excess_intervals,
      deficit_intervals = deficit_intervals
    )
    
    ## Set attributes
    attr(ret, "class") = "excess_model"
    return(ret)
  }
}

## Helper function to process identified indices into intervals
identify_intervals <- function(indices, dates) {
  if (length(indices) == 0) return(NULL)
  clusters <- cumsum(c(1, diff(indices)) > 1)
  intervals <- split(indices, clusters)
  lapply(intervals, function(x) {
    list(
      start = min(dates[x]),
      end = max(dates[x]),
      count = length(x)
    )
  })
}