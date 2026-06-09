#' Print Method for bbnp_density Objects
#'
#' @param x An object of class bbnp_density
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns the input object
#' @export
print.bbnp_density <- function(x, digits = 4, ...) {
  cat("Bias-Bounded Density Estimation\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat(sprintf("Sample size: n = %d\n", x$n))

  # Bandwidth info
  h_method <- if (is.null(x$call$h)) {
    sprintf("h = %.4f (automatic)", x$bandwidth)
  } else {
    sprintf("h = %.4f (user-specified)", x$bandwidth)
  }
  cat("Bandwidth:  ", h_method, "\n")

  cat("Kernel:     ", x$kernel, "\n\n")

  # Bias bound parameters
  cat("Bias bound parameters:\n")
  cat(sprintf("  A = %.4f, r = %.4f\n",
              x$bias_bound$est_A, x$bias_bound$est_r))
  cat(sprintf("  bias bound b1x = %.4f\n\n", x$bias_bound$b1x))

  # Estimation info
  if (!is.null(x$density)) {
    # Range estimation
    cat(sprintf("Evaluation points: %d (range: [%.4f, %.4f])\n",
                length(x$x), min(x$x), max(x$x)))
  } else if (!is.null(x$estimate)) {
    # Point estimation
    cat(sprintf("Point estimate at x = %.4f: f(x) = %.4f\n",
                x$x, x$estimate))
  }

  # Confidence level
  alpha <- x$conf_int$conf_level
  if (!is.null(alpha)) {
    cat(sprintf("Confidence level: %.0f%%\n\n",
                (1 - alpha) * 100))
  }

  cat("Use summary() for detailed statistics\n")
  cat("Use plot() to visualize results\n")

  invisible(x)
}


#' Print Method for bbnp_regression Objects
#'
#' @param x An object of class bbnp_regression
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns the input object
#' @export
print.bbnp_regression <- function(x, digits = 4, ...) {
  cat("Bias-Bounded Conditional Expectation Estimation\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat(sprintf("Sample size: n = %d\n", x$n))

  # Bandwidth info
  h_method <- if (is.null(x$call$h)) {
    sprintf("h = %.4f (automatic)", x$bandwidth)
  } else {
    sprintf("h = %.4f (user-specified)", x$bandwidth)
  }
  cat("Bandwidth:  ", h_method, "\n")

  cat("Kernel:     ", x$kernel, "\n\n")

  # Bias bound parameters
  cat("Bias bound parameters:\n")
  cat(sprintf("  A = %.4f, r = %.4f, B = %.4f\n",
              x$bias_bound$est_A, x$bias_bound$est_r, x$bias_bound$est_B))
  cat(sprintf("  bias bounds: b1x = %.4f, byx = %.4f\n\n",
              x$bias_bound$b1x, x$bias_bound$byx))

  # Estimation info
  if (!is.null(x$fitted_values)) {
    # Range estimation
    cat(sprintf("Evaluation points: %d (range: [%.4f, %.4f])\n",
                length(x$x), min(x$x), max(x$x)))
    cat(sprintf("Fitted values: E[Y|X] range [%.4f, %.4f]\n",
                min(x$fitted_values, na.rm = TRUE),
                max(x$fitted_values, na.rm = TRUE)))
  } else if (!is.null(x$estimate)) {
    # Point estimation
    cat(sprintf("Point estimate at x = %.4f: E[Y|X=x] = %.4f\n",
                x$x, x$estimate))
  }

  # Confidence level
  alpha <- x$conf_int$conf_level
  if (!is.null(alpha)) {
    cat(sprintf("Confidence level: %.0f%%\n\n",
                (1 - alpha) * 100))
  }

  cat("Use summary() for detailed statistics\n")
  cat("Use plot() to visualize results\n")
  cat("Use fitted() to extract fitted values\n")

  invisible(x)
}
