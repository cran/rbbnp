#' Plot the Fourier Transform (Deprecated)
#'
#' @description
#' **Deprecated:** This function is deprecated and will be removed in a future version.
#'
#' Use \code{plot(fit, type = "ft")} instead, where
#' \code{fit} is a \code{bbnp_density} or \code{bbnp_regression} object.
#'
#' @param X A numerical vector of sample data.
#' @param xi_interval A list containing the lower (`xi_lb`) and upper (`xi_ub`) bounds of the xi interval.
#' @param ft_plot.resol An integer representing the resolution of the plot, specifically the number of points
#'        used to represent the Fourier transform. Defaults to 500.
#'
#' @return A ggplot object representing the plot of the Fourier transform.
#'
#' @details
#' C = 1, the parameter in \eqn{O(1/n^{0.25})}, see more details in in Schennach (2020).
#'
#' @export
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Old (deprecated):
#' # plot_ft(sample_data$X, xi_interval = list(xi_lb = 1, xi_ub = 50))
#'
#' # New (recommended):
#' fit <- biasBound_density(sample_data$X, h = 0.1)
#' plot(fit, type = "ft")
#' }
plot_ft <- function(X,
                    xi_interval,
                    ft_plot.resol = 500) {

  .Deprecated(
    msg = paste0(
      "'plot_ft()' is deprecated. ",
      "Use plot(fit, type = \"ft\") instead, where fit is a bbnp_density or bbnp_regression object. ",
      "See ?biasBound_density for examples."
    )
  )

  xi_lb <- xi_interval$xi_lb
  xi_ub <- xi_interval$xi_ub

  xi_ci <- seq(xi_lb, xi_ub, length.out = log(xi_ub / xi_lb) * ft_plot.resol)

  avg_phi <- purrr::map_dbl(xi_ci, function(u) get_avg_phi(xi = u, X = X))

  log_xi = log(xi_ci)
  log_avg_phi = log(avg_phi)
  dt <- data.frame(log_xi, log_avg_phi)
  ggplot(dt, aes(x = log_xi, y = log_avg_phi)) +
    geom_line() +
    geom_vline(xintercept = c(log(xi_lb), log(xi_ub)), color = "grey") +
    annotate("text", x = max(dt$log_xi) * 0.9, y = max(dt$log_avg_phi), label = paste("N =", length(X)), size = 4)
}

#' Bias bound approach for density estimation
#'
#' Estimates the density at a given point or across a range, and provides visualization options for density,
#' bias, and confidence intervals.
#'
#' @param X A numerical vector of sample data.
#' @param x Optional. A scalar or range of points where the density is estimated. If NULL, a range is automatically generated.
#' @param h A scalar bandwidth parameter. If NULL, the bandwidth is automatically selected using the method specified in 'h_method'.
#' @param h_method Method for automatic bandwidth selection when h is NULL. Options are "cv" (cross-validation) and "silverman" (Silverman's rule of thumb). Default is "cv".
#' @param alpha Confidence level for intervals. Default is 0.05.
#' @param resol Resolution for the estimation range. Default is 100.
#' @param xi_lb Optional. Lower bound for the interval of Fourier Transform frequency xi. Used for determining the range over which A and r is estimated. If NULL, it is automatically determined based on the methods_get_xi.
#' @param xi_ub Optional. Upper bound for the interval of Fourier Transform frequency xi. Similar to xi_lb, it defines the upper range for A and r estimation. If NULL, the upper bound is determined based on the methods_get_xi.
#' @param methods_get_xi A string selecting the frequency-window rule used when xi_lb/xi_ub are NULL: "snr" (default; a signal-to-noise cutoff that selects a valid window at realistic sample sizes), "Schennach" (the data-driven rule of Schennach 2020, Theorem 2), or "Schennach_loose" (the initial, un-refined interval).
#' @param noise_floor Noise-floor form for the Schennach test: "auto" (default), "compact", or "general".
#' @param envelope_use_Y If TRUE (default), fit the regression envelope to the cross-spectrum `|phi_YX|`; if FALSE, fit it to the marginal spectrum `|phi_X|`.
#' @param integer_r If TRUE (default), clamp the fitted envelope slope up to r = 2 when it falls below the minimum smoothness assumed by Schennach (2020, Definition 2), i.e. r < 2, and refit A; this keeps the bias-bound integral finite. Slopes >= 2 are left unchanged.
#' @param ora_Ar Optional list of oracle values for A and r (for research/comparison purposes).
#' @param kernel.fun A string specifying the kernel function to be used. Options are "Schennach2004", "sinc", "normal", "epanechnikov".
#' @param if_approx_kernel Logical. If TRUE, uses approximations for the kernel function.
#' @param kernel.resol The resolution for kernel function approximation. See \code{\link{fun_approx}}.
#' @return An object of class \code{bbnp_density} with components:
#'   \item{density}{Density estimates (for range estimation)}
#'   \item{x}{Evaluation points}
#'   \item{estimate}{Point estimate (for single x)}
#'   \item{conf_int}{List containing lower, upper bounds and conf_level}
#'   \item{bias_bound}{List containing b1x, est_A, est_r, xi_interval}
#'   \item{std_error}{Standard errors}
#'   \item{call}{The function call}
#'   \item{bandwidth}{Bandwidth used}
#'   \item{n}{Sample size}
#'   \item{kernel}{Kernel type}
#'   \item{data}{Original data}
#'
#'   Use \code{plot()}, \code{summary()}, \code{coef()}, and \code{confint()} methods to work with the result.
#' @export
#' @examples
#' \donttest{
#' # Example 1: Point estimation at x = 1
#' X <- rnorm(100)
#' fit <- biasBound_density(X = X, x = 1, h = 0.09)
#' print(fit)
#' coef(fit)
#'
#' # Example 2: Range estimation with automatic bandwidth
#' fit2 <- biasBound_density(X = X, h = NULL, h_method = "cv")
#' plot(fit2)           # Density plot
#' plot(fit2, type = "ft")  # Fourier transform plot
#' summary(fit2)
#' }
biasBound_density <- function(X, x = NULL, h = NULL, h_method = "cv", alpha = 0.05, resol = 100,
                              xi_lb = NULL, xi_ub = NULL, methods_get_xi = "snr",
                              noise_floor = "auto", envelope_use_Y = TRUE, integer_r = TRUE,
                              ora_Ar = NULL,
                              kernel.fun = "Schennach2004", if_approx_kernel = TRUE, kernel.resol = 1000) {

  # Capture the function call
  call <- match.call()

  # Create a configuration object that contains all settings and pre-computed values
  config <- create_biasBound_config(
    X = X,
    h = h,
    h_method = h_method,
    alpha = alpha,
    resol = resol,
    xi_lb = xi_lb,
    xi_ub = xi_ub,
    methods_get_xi = methods_get_xi,
    noise_floor = noise_floor,
    envelope_use_Y = envelope_use_Y,
    integer_r = integer_r,
    kernel.fun = kernel.fun,
    if_approx_kernel = if_approx_kernel,
    kernel.resol = kernel.resol
  )

  # Extract necessary values from config
  inf_k <- config$kernel_functions$kernel
  xi_interval <- config$xi_interval
  est_Ar <- config$est_Ar
  b1x <- config$b1x
  h <- config$h  # Get the possibly auto-selected bandwidth

  # Compute estimates based on whether x is specified (point vs range estimation)
  if (is.null(x)) {
    # Range estimation
    x_range <- seq(min(X) - sd(X)*0.5, max(X) + sd(X) * 0.5, length.out = resol)
    f1x <- purrr::map(x_range, get_avg_f1x, X = X, h = h, inf_k = inf_k) %>% unlist()
    f1x <- pmax(f1x, 0)
    sigma <- purrr::map(x_range, get_sigma, X = X, h = h, inf_k = inf_k) %>% unlist()

    # Compute confidence intervals
    lb <- pmax(f1x - sigma * qnorm(1 - alpha / 2) - b1x, 0)
    ub <- pmax(f1x + sigma * qnorm(1 - alpha / 2) + b1x, 0)

    # Create S3 object for range estimation
    result <- new_bbnp_density(
      density = f1x,
      x = x_range,
      estimate = NULL,
      conf_int = list(
        lower = lb,
        upper = ub,
        conf_level = alpha
      ),
      bias_bound = list(
        b1x = b1x,
        est_A = est_Ar[1],
        est_r = est_Ar[2],
        xi_interval = xi_interval
      ),
      std_error = sigma,
      call = call,
      bandwidth = h,
      n = length(X),
      kernel = kernel.fun,
      data = list(X = X),
      internals = list(
        config = config,
        kernel_functions = config$kernel_functions
      )
    )

  } else {
    # Point estimation for specific x value
    f1x <- get_avg_f1x(X, x, h, inf_k = inf_k)
    f1x <- max(f1x, 0)
    sigma <- get_sigma(X, x, h, inf_k = inf_k)

    # Compute confidence interval
    lb <- max(c(f1x - sigma * qnorm(1 - alpha / 2) - b1x, 0))
    ub <- max(c(f1x + sigma * qnorm(1 - alpha / 2) + b1x, 0))

    # Create S3 object for point estimation
    result <- new_bbnp_density(
      density = NULL,
      x = x,
      estimate = f1x,
      conf_int = list(
        lower = lb,
        upper = ub,
        conf_level = alpha
      ),
      bias_bound = list(
        b1x = b1x,
        est_A = est_Ar[1],
        est_r = est_Ar[2],
        xi_interval = xi_interval
      ),
      std_error = sigma,
      call = call,
      bandwidth = h,
      n = length(X),
      kernel = kernel.fun,
      data = list(X = X),
      internals = list(
        config = config,
        kernel_functions = config$kernel_functions
      )
    )
  }

  return(result)
}

#' Bias bound approach for conditional expectation estimation
#'
#' Estimates the density at a given point or across a range, and provides visualization options for density,
#' bias, and confidence intervals.
#'
#' @param Y A numerical vector of sample data.
#' @param X A numerical vector of sample data.
#' @param x Optional. A scalar or range of points where the density is estimated. If NULL, a range is automatically generated.
#' @param h A scalar bandwidth parameter. If NULL, the bandwidth is automatically selected using the method specified in 'h_method'.
#' @param h_method Method for automatic bandwidth selection when h is NULL. Options are "cv" (cross-validation) and "silverman" (Silverman's rule of thumb). Default is "cv".
#' @param alpha Confidence level for intervals. Default is 0.05.
#' @param est_Ar Optional list of estimates for A and r. If NULL, they are computed using `get_est_Ar()`.
#' @param resol Resolution for the estimation range. Default is 100.
#' @param xi_lb Optional. Lower bound for the interval of Fourier Transform frequency xi. Used for determining the range over which A and r is estimated. If NULL, it is automatically determined based on the methods_get_xi.
#' @param xi_ub Optional. Upper bound for the interval of Fourier Transform frequency xi. Similar to xi_lb, it defines the upper range for A and r estimation. If NULL, the upper bound is determined based on the methods_get_xi.
#' @param methods_get_xi A string selecting the frequency-window rule used when xi_lb/xi_ub are NULL: "snr" (default; a signal-to-noise cutoff that selects a valid window at realistic sample sizes), "Schennach" (the data-driven rule of Schennach 2020, Theorem 2), or "Schennach_loose" (the initial, un-refined interval).
#' @param noise_floor Noise-floor form for the Schennach test: "auto" (default), "compact", or "general".
#' @param envelope_use_Y If TRUE (default), fit the regression envelope to the cross-spectrum `|phi_YX|`; if FALSE, fit it to the marginal spectrum `|phi_X|`.
#' @param integer_r If TRUE (default), clamp the fitted envelope slope up to r = 2 when it falls below the minimum smoothness assumed by Schennach (2020, Definition 2), i.e. r < 2, and refit A; this keeps the bias-bound integral finite. Slopes >= 2 are left unchanged.
#' @param ora_Ar Optional list of oracle values for A and r (for research/comparison purposes).
#' @param kernel.fun A string specifying the kernel function to be used. Options are "Schennach2004", "sinc", "normal", "epanechnikov".
#' @param if_approx_kernel Logical. If TRUE, uses approximations for the kernel function.
#' @param kernel.resol The resolution for kernel function approximation. See \code{\link{fun_approx}}.
#' @return An object of class \code{bbnp_regression} with components:
#'   \item{fitted_values}{\eqn{E[Y|X=x]} estimates (for range estimation)}
#'   \item{x}{Evaluation points}
#'   \item{estimate}{Point estimate (for single x)}
#'   \item{conf_int}{List containing lower, upper bounds and conf_level. Note that
#'     the confidence interval can be unbounded (i.e., contain \code{-Inf} or \code{Inf})
#'     in regions where the estimated marginal density \eqn{\hat f(x)} is very close to zero,
#'     because the estimator is formed as a ratio involving \eqn{1/\hat f(x)}.}
#'   \item{bias_bound}{List containing b1x, byx, est_A, est_r, est_B, xi_interval}
#'   \item{std_error}{Standard errors}
#'   \item{marginal_density}{f(x) estimates}
#'   \item{joint_density}{f_YX estimates}
#'   \item{call}{The function call}
#'   \item{bandwidth}{Bandwidth used}
#'   \item{n}{Sample size}
#'   \item{kernel}{Kernel type}
#'   \item{data}{Original data (X, Y)}
#'
#'   Use \code{plot()}, \code{summary()}, \code{coef()}, \code{fitted()}, and \code{confint()} methods to work with the result.
#' @export
#' @examples
#' \donttest{
#' # Example 1: Point estimation at x = 1
#' X <- rnorm(100)
#' Y <- X^2 + rnorm(100)
#' fit <- biasBound_condExpectation(Y = Y, X = X, x = 1, h = 0.09)
#' print(fit)
#' fitted(fit)
#'
#' # Example 2: Range estimation with plots
#' fit2 <- biasBound_condExpectation(Y = Y, X = X, h = NULL, h_method = "cv")
#' plot(fit2)              # Regression plot
#' plot(fit2, type = "ft") # Fourier transform plot
#' summary(fit2)
#' }
biasBound_condExpectation <- function(Y, X, x = NULL, h = NULL, h_method = "cv", alpha = 0.05, est_Ar = NULL, resol = 100,
                                      xi_lb = NULL, xi_ub = NULL, methods_get_xi = "snr",
                                      noise_floor = "auto", envelope_use_Y = TRUE, integer_r = TRUE,
                                      ora_Ar = NULL,
                                      kernel.fun = "Schennach2004", if_approx_kernel = TRUE, kernel.resol = 1000) {
  # Capture the function call
  call <- match.call()

  # regularization of Y and X data structure
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length!")
  }

  # Create a configuration object that contains all settings and pre-computed values
  config <- create_biasBound_config(
    X = X,
    Y = Y,
    h = h,
    h_method = h_method,
    alpha = alpha,
    resol = resol,
    xi_lb = xi_lb,
    xi_ub = xi_ub,
    methods_get_xi = methods_get_xi,
    noise_floor = noise_floor,
    envelope_use_Y = envelope_use_Y,
    integer_r = integer_r,
    kernel.fun = kernel.fun,
    if_approx_kernel = if_approx_kernel,
    kernel.resol = kernel.resol
  )

  # Extract necessary values from config
  inf_k <- config$kernel_functions$kernel
  xi_interval <- config$xi_interval
  est_Ar <- config$est_Ar
  est_Ar_1x <- config$est_Ar_1x
  b1x <- config$b1x
  est_B <- config$est_B
  byx <- config$byx
  h <- config$h  # Get the possibly auto-selected bandwidth

  # Compute estimates based on whether x is specified (point vs range estimation)
  if (is.null(x)) {
    # Range estimation
    x_range <- seq(min(X), max(X), length.out = resol)
    f1x <- purrr::map(x_range, get_avg_f1x, X = X, h = h, inf_k = inf_k) %>% unlist()
    fyx <- purrr::map(x_range, get_avg_fyx, Y = Y, X = X, h = h, inf_k = inf_k) %>% unlist()
    sigma <- purrr::map(x_range, get_sigma, X = X, h = h, inf_k = inf_k) %>% unlist()
    sigma_yx <- purrr::map(x_range, get_sigma_yx, Y = Y, X = X, h = h, inf_k = inf_k) %>% unlist()

    # Compute confidence intervals
    lb <- (fyx - byx) / pmax(f1x + sign(fyx - byx) * b1x, 0) - sigma_yx * qnorm(1 - alpha / 2)
    ub <- (fyx + byx) / pmax(f1x - sign(fyx + byx) * b1x, 0) + sigma_yx * qnorm(1 - alpha / 2)

    conditional_mean_yx <- fyx / f1x

    # Create S3 object for range estimation
    result <- new_bbnp_regression(
      fitted_values = conditional_mean_yx,
      x = x_range,
      estimate = NULL,
      conf_int = list(
        lower = lb,
        upper = ub,
        conf_level = alpha
      ),
      bias_bound = list(
        b1x = b1x,
        byx = byx,

        # Envelope parameters used for the joint object (Y;X)
        est_A = est_Ar[1],
        est_r = est_Ar[2],

        # Envelope parameters used for the marginal object (1;X) that underlies b1x
        est_A_1x = est_Ar_1x[1],
        est_r_1x = est_Ar_1x[2],

        est_B = est_B,
        xi_interval = xi_interval
      ),
      std_error = sigma_yx,
      marginal_density = f1x,
      joint_density = fyx,
      call = call,
      bandwidth = h,
      n = length(X),
      kernel = kernel.fun,
      data = list(X = X, Y = Y),
      internals = list(
        config = config,
        kernel_functions = config$kernel_functions
      )
    )

  } else {
    # Point estimation for specific x value
    f1x <- get_avg_f1x(X, x, h, inf_k = inf_k)
    fyx <- get_avg_fyx(Y = Y, X = X, x = x, h = h, inf_k = inf_k)
    sigma <- get_sigma(X, x, h, inf_k = inf_k)
    sigma_yx <- get_sigma_yx(Y = Y, X = X, x = x, h = h, inf_k = inf_k)
    conditional_mean_yx <- fyx / f1x

    # Compute confidence interval
    lb <- (fyx - byx) / max(c(f1x + sign(fyx - byx) * b1x, 0)) - sigma_yx * qnorm(1 - alpha / 2)
    ub <- (fyx + byx) / max(c(f1x - sign(fyx + byx) * b1x, 0)) + sigma_yx * qnorm(1 - alpha / 2)

    # Create S3 object for point estimation
    result <- new_bbnp_regression(
      fitted_values = NULL,
      x = x,
      estimate = conditional_mean_yx,
      conf_int = list(
        lower = lb,
        upper = ub,
        conf_level = alpha
      ),
      bias_bound = list(
        b1x = b1x,
        byx = byx,

        # Envelope parameters used for the joint object (Y;X)
        est_A = est_Ar[1],
        est_r = est_Ar[2],

        # Envelope parameters used for the marginal object (1;X) that underlies b1x
        est_A_1x = est_Ar_1x[1],
        est_r_1x = est_Ar_1x[2],

        est_B = est_B,
        xi_interval = xi_interval
      ),
      std_error = sigma_yx,
      marginal_density = f1x,
      joint_density = fyx,
      call = call,
      bandwidth = h,
      n = length(X),
      kernel = kernel.fun,
      data = list(X = X, Y = Y),
      internals = list(
        config = config,
        kernel_functions = config$kernel_functions
      )
    )
  }

  return(result)
}
