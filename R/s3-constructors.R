#' S3 Constructor for bbnp_density Objects
#'
#' Creates a bbnp_density S3 object with validation
#'
#' @param density Numeric vector of density estimates (for range estimation)
#' @param x Numeric vector of evaluation points
#' @param estimate Numeric scalar point estimate (for point estimation)
#' @param conf_int List containing lower, upper, and conf_level
#' @param bias_bound List containing b1x, est_A, est_r, xi_interval
#' @param std_error Numeric vector of standard errors
#' @param call The original function call
#' @param bandwidth Numeric scalar bandwidth used
#' @param n Integer sample size
#' @param kernel Character string kernel type
#' @param data List containing original data (X)
#' @param internals List containing internal objects (config, kernel_functions)
#'
#' @return An object of class c("bbnp_density", "bbnp")
#' @keywords internal
new_bbnp_density <- function(density = NULL,
                              x = NULL,
                              estimate = NULL,
                              conf_int = list(),
                              bias_bound = list(),
                              std_error = NULL,
                              call = NULL,
                              bandwidth = NULL,
                              n = NULL,
                              kernel = NULL,
                              data = list(),
                              internals = list()) {

  # Validation
  if (!is.null(density) && !is.numeric(density)) {
    stop("density must be numeric", call. = FALSE)
  }
  if (!is.null(x) && !is.numeric(x)) {
    stop("x must be numeric", call. = FALSE)
  }
  if (!is.null(estimate) && (!is.numeric(estimate) || length(estimate) != 1)) {
    stop("estimate must be a numeric scalar", call. = FALSE)
  }
  if (!is.null(bandwidth) && (!is.numeric(bandwidth) || length(bandwidth) != 1 || bandwidth <= 0)) {
    stop("bandwidth must be a positive numeric scalar", call. = FALSE)
  }
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1 || n <= 0)) {
    stop("n must be a positive integer", call. = FALSE)
  }

  # Build structure
  structure(
    list(
      density = density,
      x = x,
      estimate = estimate,
      conf_int = conf_int,
      bias_bound = bias_bound,
      std_error = std_error,
      call = call,
      bandwidth = bandwidth,
      n = n,
      kernel = kernel,
      data = data,
      internals = internals,
      .internals = internals
    ),
    class = c("bbnp_density", "bbnp")
  )
}


#' S3 Constructor for bbnp_regression Objects
#'
#' Creates a bbnp_regression S3 object with validation
#'
#' @param fitted_values Numeric vector of \eqn{E[Y|X=x]} estimates
#' @param x Numeric vector of evaluation points
#' @param estimate Numeric scalar point estimate (for point estimation)
#' @param conf_int List containing lower, upper, and conf_level
#' @param bias_bound List containing b1x, byx, est_A, est_r, est_B, xi_interval
#' @param std_error Numeric vector of standard errors
#' @param marginal_density Numeric vector of f(x) estimates
#' @param joint_density Numeric vector of f_YX estimates
#' @param call The original function call
#' @param bandwidth Numeric scalar bandwidth used
#' @param n Integer sample size
#' @param kernel Character string kernel type
#' @param data List containing original data (X, Y)
#' @param internals List containing internal objects (config, kernel_functions)
#'
#' @return An object of class c("bbnp_regression", "bbnp")
#' @keywords internal
new_bbnp_regression <- function(fitted_values = NULL,
                                 x = NULL,
                                 estimate = NULL,
                                 conf_int = list(),
                                 bias_bound = list(),
                                 std_error = NULL,
                                 marginal_density = NULL,
                                 joint_density = NULL,
                                 call = NULL,
                                 bandwidth = NULL,
                                 n = NULL,
                                 kernel = NULL,
                                 data = list(),
                                 internals = list()) {

  # Validation
  if (!is.null(fitted_values) && !is.numeric(fitted_values)) {
    stop("fitted_values must be numeric", call. = FALSE)
  }
  if (!is.null(x) && !is.numeric(x)) {
    stop("x must be numeric", call. = FALSE)
  }
  if (!is.null(estimate) && (!is.numeric(estimate) || length(estimate) != 1)) {
    stop("estimate must be a numeric scalar", call. = FALSE)
  }
  if (!is.null(bandwidth) && (!is.numeric(bandwidth) || length(bandwidth) != 1 || bandwidth <= 0)) {
    stop("bandwidth must be a positive numeric scalar", call. = FALSE)
  }
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1 || n <= 0)) {
    stop("n must be a positive integer", call. = FALSE)
  }
  if (!is.null(marginal_density) && !is.numeric(marginal_density)) {
    stop("marginal_density must be numeric", call. = FALSE)
  }
  if (!is.null(joint_density) && !is.numeric(joint_density)) {
    stop("joint_density must be numeric", call. = FALSE)
  }

  # Build structure
  structure(
    list(
      fitted_values = fitted_values,
      x = x,
      estimate = estimate,
      conf_int = conf_int,
      bias_bound = bias_bound,
      std_error = std_error,
      marginal_density = marginal_density,
      joint_density = joint_density,
      call = call,
      bandwidth = bandwidth,
      n = n,
      kernel = kernel,
      data = data,
      internals = internals,
      .internals = internals
    ),
    class = c("bbnp_regression", "bbnp")
  )
}
