#' Summary Method for bbnp_density Objects
#'
#' @param object An object of class bbnp_density
#' @param ... Additional arguments (unused)
#'
#' @return An object of class summary.bbnp_density
#' @export
summary.bbnp_density <- function(object, ...) {
  structure(
    list(
      call = object$call,
      n = object$n,
      bandwidth = object$bandwidth,
      kernel = object$kernel,
      bias_bound = object$bias_bound,
      conf_int = object$conf_int,
      is_point_est = !is.null(object$estimate),
      estimate = object$estimate,
      x = object$x,
      density_summary = if (!is.null(object$density)) {
        c(
          min = min(object$density, na.rm = TRUE),
          Q1 = quantile(object$density, 0.25, na.rm = TRUE),
          median = median(object$density, na.rm = TRUE),
          mean = mean(object$density, na.rm = TRUE),
          Q3 = quantile(object$density, 0.75, na.rm = TRUE),
          max = max(object$density, na.rm = TRUE)
        )
      } else NULL,
      std_error_summary = if (!is.null(object$std_error)) {
        if (length(object$std_error) > 1) {
          c(
            min = min(object$std_error, na.rm = TRUE),
            mean = mean(object$std_error, na.rm = TRUE),
            max = max(object$std_error, na.rm = TRUE)
          )
        } else {
          object$std_error
        }
      } else NULL
    ),
    class = "summary.bbnp_density"
  )
}


#' Print Method for summary.bbnp_density Objects
#'
#' @param x An object of class summary.bbnp_density
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns the input object
#' @keywords internal
#' @export
print.summary.bbnp_density <- function(x, digits = 4, ...) {
  cat("Summary: Bias-Bounded Density Estimation\n")
  cat(rep("=", 60), "\n\n", sep = "")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Sample Information:\n")
  cat(sprintf("  Sample size (n):  %d\n", x$n))
  cat(sprintf("  Bandwidth (h):    %.4f\n", x$bandwidth))
  cat(sprintf("  Kernel function:  %s\n\n", x$kernel))

  cat("Bias Bound Parameters:\n")
  cat(sprintf("  A (amplitude):    %.4f\n", x$bias_bound$est_A))
  cat(sprintf("  r (decay rate):   %.4f\n", x$bias_bound$est_r))
  cat(sprintf("  b1x (bias bound): %.4f\n", x$bias_bound$b1x))
  if (!is.null(x$bias_bound$xi_interval)) {
    cat(sprintf("  Xi interval:      [%.4f, %.4f]\n",
                x$bias_bound$xi_interval$xi_lb,
                x$bias_bound$xi_interval$xi_ub))
  }
  cat("\n")

  if (x$is_point_est) {
    cat("Point Estimation:\n")
    cat(sprintf("  x = %.4f\n", x$x))
    cat(sprintf("  f(x) = %.4f\n", x$estimate))
    if (!is.null(x$conf_int$lower)) {
      cat(sprintf("  %d%% CI: [%.4f, %.4f]\n",
                  (1 - x$conf_int$conf_level) * 100,
                  x$conf_int$lower,
                  x$conf_int$upper))
    }
  } else {
    cat("Range Estimation:\n")
    cat("  Density estimates:\n")
    print(round(x$density_summary, digits))
    cat("\n")
    if (!is.null(x$std_error_summary)) {
      cat("  Standard errors:\n")
      print(round(x$std_error_summary, digits))
      cat("\n")
    }
  }

  invisible(x)
}


#' Summary Method for bbnp_regression Objects
#'
#' @param object An object of class bbnp_regression
#' @param ... Additional arguments (unused)
#'
#' @return An object of class summary.bbnp_regression
#' @export
summary.bbnp_regression <- function(object, ...) {
  structure(
    list(
      call = object$call,
      n = object$n,
      bandwidth = object$bandwidth,
      kernel = object$kernel,
      bias_bound = object$bias_bound,
      conf_int = object$conf_int,
      is_point_est = !is.null(object$estimate),
      estimate = object$estimate,
      x = object$x,
      fitted_summary = if (!is.null(object$fitted_values)) {
        c(
          min = min(object$fitted_values, na.rm = TRUE),
          Q1 = quantile(object$fitted_values, 0.25, na.rm = TRUE),
          median = median(object$fitted_values, na.rm = TRUE),
          mean = mean(object$fitted_values, na.rm = TRUE),
          Q3 = quantile(object$fitted_values, 0.75, na.rm = TRUE),
          max = max(object$fitted_values, na.rm = TRUE)
        )
      } else NULL,
      marginal_density_summary = if (!is.null(object$marginal_density) && length(object$marginal_density) > 1) {
        c(
          min = min(object$marginal_density, na.rm = TRUE),
          mean = mean(object$marginal_density, na.rm = TRUE),
          max = max(object$marginal_density, na.rm = TRUE)
        )
      } else object$marginal_density,
      std_error_summary = if (!is.null(object$std_error)) {
        if (length(object$std_error) > 1) {
          c(
            min = min(object$std_error, na.rm = TRUE),
            mean = mean(object$std_error, na.rm = TRUE),
            max = max(object$std_error, na.rm = TRUE)
          )
        } else {
          object$std_error
        }
      } else NULL
    ),
    class = "summary.bbnp_regression"
  )
}


#' Print Method for summary.bbnp_regression Objects
#'
#' @param x An object of class summary.bbnp_regression
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns the input object
#' @keywords internal
#' @export
print.summary.bbnp_regression <- function(x, digits = 4, ...) {
  cat("Summary: Bias-Bounded Conditional Expectation Estimation\n")
  cat(rep("=", 60), "\n\n", sep = "")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Sample Information:\n")
  cat(sprintf("  Sample size (n):  %d\n", x$n))
  cat(sprintf("  Bandwidth (h):    %.4f\n", x$bandwidth))
  cat(sprintf("  Kernel function:  %s\n\n", x$kernel))

  cat("Bias Bound Parameters:\n")
  cat(sprintf("  A (amplitude):    %.4f\n", x$bias_bound$est_A))
  cat(sprintf("  r (decay rate):   %.4f\n", x$bias_bound$est_r))
  cat(sprintf("  B (Y bound):      %.4f\n", x$bias_bound$est_B))
  cat(sprintf("  b1x (bias f(x)):  %.4f\n", x$bias_bound$b1x))
  cat(sprintf("  byx (bias f_YX):  %.4f\n", x$bias_bound$byx))
  if (!is.null(x$bias_bound$xi_interval)) {
    cat(sprintf("  Xi interval:      [%.4f, %.4f]\n",
                x$bias_bound$xi_interval$xi_lb,
                x$bias_bound$xi_interval$xi_ub))
  }
  cat("\n")

  if (x$is_point_est) {
    cat("Point Estimation:\n")
    cat(sprintf("  x = %.4f\n", x$x))
    cat(sprintf("  E[Y|X=x] = %.4f\n", x$estimate))
    if (!is.null(x$conf_int$lower)) {
      cat(sprintf("  %d%% CI: [%.4f, %.4f]\n",
                  (1 - x$conf_int$conf_level) * 100,
                  x$conf_int$lower,
                  x$conf_int$upper))
    }
  } else {
    cat("Range Estimation:\n")
    cat("  Fitted values (E[Y|X]):\n")
    print(round(x$fitted_summary, digits))
    cat("\n")
    if (!is.null(x$marginal_density_summary) && length(x$marginal_density_summary) > 1) {
      cat("  Marginal density f(x):\n")
      print(round(x$marginal_density_summary, digits))
      cat("\n")
    }
    if (!is.null(x$std_error_summary)) {
      cat("  Standard errors:\n")
      print(round(x$std_error_summary, digits))
      cat("\n")
    }
  }

  invisible(x)
}
