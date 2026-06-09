#' Extract Coefficients from bbnp_density Object
#'
#' Extracts the estimated bias bound parameters and bandwidth
#'
#' @param object An object of class bbnp_density
#' @param ... Additional arguments (unused)
#'
#' @return Named numeric vector with A (amplitude), r (decay rate), and h (bandwidth)
#' @export
#'
#' @examples
#' \donttest{
#' X <- rnorm(100)
#' fit <- biasBound_density(X, h = 0.1)
#' coef(fit)
#' }
coef.bbnp_density <- function(object, ...) {
  c(
    A = unname(object$bias_bound$est_A),
    r = unname(object$bias_bound$est_r),
    h = unname(object$bandwidth)
  )
}


#' Extract Coefficients from bbnp_regression Object
#'
#' Extracts the estimated bias bound parameters and bandwidth
#'
#' @param object An object of class bbnp_regression
#' @param ... Additional arguments (unused)
#'
#' @return Named numeric vector with A (amplitude), r (decay rate), B (Y bound), and h (bandwidth)
#' @export
#'
#' @examples
#' \donttest{
#' X <- rnorm(100)
#' Y <- X^2 + rnorm(100)
#' fit <- biasBound_condExpectation(Y, X, h = 0.1)
#' coef(fit)
#' }
coef.bbnp_regression <- function(object, ...) {
  c(
    A = unname(object$bias_bound$est_A),
    r = unname(object$bias_bound$est_r),
    B = unname(object$bias_bound$est_B),
    h = unname(object$bandwidth)
  )
}


#' Extract Fitted Values from bbnp_regression Object
#'
#' Extracts the fitted conditional expectation values \eqn{E[Y|X=x]}
#'
#' @param object An object of class bbnp_regression
#' @param ... Additional arguments (unused)
#'
#' @return Numeric vector of fitted values for range estimation,
#'   or a single numeric value for point estimation
#' @export
#'
#' @examples
#' \donttest{
#' X <- rnorm(100)
#' Y <- X^2 + rnorm(100)
#' fit <- biasBound_condExpectation(Y, X, h = 0.1)
#' fitted(fit)
#' }
fitted.bbnp_regression <- function(object, ...) {
  if (!is.null(object$fitted_values)) {
    return(object$fitted_values)
  } else if (!is.null(object$estimate)) {
    return(object$estimate)
  } else {
    stop("No fitted values available in object", call. = FALSE)
  }
}


#' Extract Confidence Intervals from bbnp_density Object
#'
#' Extracts confidence intervals for density estimates
#'
#' @param object An object of class bbnp_density
#' @param parm Not used (included for S3 generic compatibility)
#' @param level Confidence level (default: 0.95). Note: this parameter is not used
#'   as the confidence level is fixed at object creation time.
#' @param ... Additional arguments (unused)
#'
#' @return For range estimation: a matrix with columns "lower" and "upper"
#'   For point estimation: a named vector with elements "lower" and "upper"
#' @export
#'
#' @examples
#' \donttest{
#' X <- rnorm(100)
#' fit <- biasBound_density(X, h = 0.1)
#' confint(fit)
#' }
confint.bbnp_density <- function(object, parm = NULL, level = 0.95, ...) {
  # Note: level parameter is ignored as CI was computed at estimation time
  if (!is.null(level) && level != 0.95) {
    warning("The 'level' parameter is ignored. Confidence intervals are computed ",
            "at the alpha level specified during estimation (alpha = ",
            object$conf_int$conf_level, ").",
            call. = FALSE)
  }

  lower <- object$conf_int$lower
  upper <- object$conf_int$upper

  if (is.null(lower) || is.null(upper)) {
    stop("No confidence intervals available in object", call. = FALSE)
  }

  # For range estimation, return matrix
  if (length(lower) > 1) {
    return(cbind(lower = lower, upper = upper))
  } else {
    # For point estimation, return named vector
    return(c(lower = lower, upper = upper))
  }
}


#' Extract Confidence Intervals from bbnp_regression Object
#'
#' Extracts confidence intervals for conditional expectation estimates
#'
#' @param object An object of class bbnp_regression
#' @param parm Not used (included for S3 generic compatibility)
#' @param level Confidence level (default: 0.95). Note: this parameter is not used
#'   as the confidence level is fixed at object creation time.
#' @param ... Additional arguments (unused)
#'
#' @return For range estimation: a matrix with columns "lower" and "upper"
#'   For point estimation: a named vector with elements "lower" and "upper"
#' @export
#'
#' @examples
#' \donttest{
#' X <- rnorm(100)
#' Y <- X^2 + rnorm(100)
#' fit <- biasBound_condExpectation(Y, X, h = 0.1)
#' confint(fit)
#' }
confint.bbnp_regression <- function(object, parm = NULL, level = 0.95, ...) {
  # Note: level parameter is ignored as CI was computed at estimation time
  if (!is.null(level) && level != 0.95) {
    warning("The 'level' parameter is ignored. Confidence intervals are computed ",
            "at the alpha level specified during estimation (alpha = ",
            object$conf_int$conf_level, ").",
            call. = FALSE)
  }

  lower <- object$conf_int$lower
  upper <- object$conf_int$upper

  if (is.null(lower) || is.null(upper)) {
    stop("No confidence intervals available in object", call. = FALSE)
  }

  # For range estimation, return matrix
  if (length(lower) > 1) {
    return(cbind(lower = lower, upper = upper))
  } else {
    # For point estimation, return named vector
    return(c(lower = lower, upper = upper))
  }
}
