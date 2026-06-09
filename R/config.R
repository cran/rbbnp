#' Create a configuration object for bias bound estimations
#'
#' @param X A numerical vector of sample data.
#' @param Y Optional. A numerical vector of sample data for conditional expectation.
#' @param h A scalar bandwidth parameter. If NULL, the bandwidth is automatically selected using the method specified in 'h_method'.
#' @param h_method Method for automatic bandwidth selection when h is NULL. Options are "cv" (cross-validation) and "silverman" (Silverman's rule of thumb). Default is "cv".
#' @param use_fft Ignored. Maintained for backward compatibility.
#' @param alpha Confidence level for intervals.
#' @param resol Resolution for the estimation range.
#' @param xi_lb Lower bound for the interval of Fourier Transform frequency.
#' @param xi_ub Upper bound for the interval of Fourier Transform frequency.
#' @param methods_get_xi Method to determine the xi interval: "snr" (default), "Schennach",
#'   or "Schennach_loose".
#' @param noise_floor Noise-floor form for the Schennach test: "auto" (default), "compact",
#'   or "general".
#' @param envelope_use_Y If TRUE (default), fit the regression envelope to the cross-spectrum
#'   `|phi_YX|`; if FALSE, fit it to the marginal spectrum `|phi_X|`.
#' @param integer_r If TRUE (default), when the fitted envelope slope falls below the minimum
#'   smoothness assumed by Schennach (2020, Definition 2), i.e. r < 2, clamp it up to r = 2 and
#'   refit A; this keeps the bias-bound integral finite. Slopes >= 2 are left unchanged.
#' @param kernel.fun Kernel function to be used. Options include "normal", "epanechnikov", "Schennach2004", and "sinc".
#' @param if_approx_kernel Use approximations for the kernel function.
#' @param kernel.resol Resolution for kernel approximation.
#'
#' @return A configuration object (list) with all parameters
#'
#' @export
create_biasBound_config <- function(X,
                                   Y = NULL,
                                   h = NULL,
                                   h_method = "cv",
                                   use_fft = TRUE,
                                   alpha = 0.05,
                                   resol = 100,
                                   xi_lb = NULL,
                                   xi_ub = NULL,
                                   methods_get_xi = "snr",
                                   noise_floor = "auto",
                                   envelope_use_Y = TRUE,
                                   integer_r = TRUE,
                                   kernel.fun = "Schennach2004",
                                   if_approx_kernel = TRUE,
                                   kernel.resol = 1000) {

  # Create kernel functions
  kernel_functions <- create_kernel_functions(
    kernel.fun = kernel.fun,
    if_approx_kernel = if_approx_kernel,
    kernel.resol = kernel.resol
  )

  # Select bandwidth if not provided
  if (is.null(h)) {
    # Use appropriate bandwidth selection method
    h <- select_bandwidth(X, Y, method = h_method, kernel.fun = kernel.fun)
  }

  # Get xi interval if not specified
  if (is.null(xi_lb) || is.null(xi_ub)) {
    # For conditional expectation settings, the Schennach test depends on v_Y.
    # Use Y when available; density case defaults to Y=1.
    Y_for_xi <- if (is.null(Y)) 1 else Y
    xi_interval <- get_xi_interval(Y = Y_for_xi, X = X, methods = methods_get_xi, noise_floor = noise_floor)
    xi_lb <- xi_interval$xi_lb
    xi_ub <- xi_interval$xi_ub
  } else {
    xi_interval <- list(xi_lb = xi_lb, xi_ub = xi_ub)
  }

  # A/r for the marginal density f_{1;X} (always uses Y=1)
  est_Ar_1x <- get_est_Ar(Y = 1, X = X, xi_interval = xi_interval, integer_r = integer_r)

  # b_{1;X} uses the marginal (Y=1) envelope
  b1x <- get_est_b1x(
    X = X,
    h = h,
    est_Ar = est_Ar_1x,
    inf_k_ft = kernel_functions$kernel_ft
  )

  # If Y is provided, estimate A/r for the joint object f_{Y;X} using phi_{Y;X}.
  # Otherwise, reuse the marginal estimates.
  if (!is.null(Y)) {
    est_Ar_yx <- get_est_Ar(Y = Y, X = X, xi_interval = xi_interval, use_Y = envelope_use_Y, integer_r = integer_r)
    est_B <- get_est_B(Y = Y)
    byx <- get_est_byx(
      Y = Y,
      X = X,
      h = h,
      est_Ar = est_Ar_yx,
      est_B = est_B,
      inf_k_ft = kernel_functions$kernel_ft
    )
    est_Ar <- est_Ar_yx
  } else {
    est_Ar_yx <- NULL
    est_B <- NULL
    byx <- NULL
    est_Ar <- est_Ar_1x
  }

  # Create and return config object
  config <- list(
    X = X,
    Y = Y,
    h = h,
    h_method = h_method,
    alpha = alpha,
    resol = resol,
    xi_lb = xi_lb,
    xi_ub = xi_ub,
    methods_get_xi = methods_get_xi,
    kernel.fun = kernel.fun,
    if_approx_kernel = if_approx_kernel,
    kernel.resol = kernel.resol,
    kernel_functions = kernel_functions,
    xi_interval = xi_interval,

    # Envelope parameters
    # - est_Ar: the primary (A,r) used by the top-level estimator
    #   * density: (A,r) for phi_{1;X}
    #   * regression: (A,r) for phi_{Y;X}
    # - est_Ar_1x: always (A,r) for phi_{1;X}
    # - est_Ar_yx: (A,r) for phi_{Y;X} (NULL when Y is NULL)
    est_Ar = est_Ar,
    est_Ar_1x = est_Ar_1x,
    est_Ar_yx = est_Ar_yx,

    b1x = b1x,
    est_B = est_B,
    byx = byx
  )

  return(config)
}
