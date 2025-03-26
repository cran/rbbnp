#' Create a configuration object for bias bound estimations
#'
#' @param X A numerical vector of sample data.
#' @param Y Optional. A numerical vector of sample data for conditional expectation.
#' @param h A scalar bandwidth parameter.
#' @param alpha Confidence level for intervals.
#' @param resol Resolution for the estimation range.
#' @param xi_lb Lower bound for the interval of Fourier Transform frequency.
#' @param xi_ub Upper bound for the interval of Fourier Transform frequency.
#' @param methods_get_xi Method to determine xi interval.
#' @param kernel.fun Kernel function to be used.
#' @param if_approx_kernel Use approximations for the kernel function.
#' @param kernel.resol Resolution for kernel approximation.
#'
#' @return A configuration object (list) with all parameters
#'
#' @export
create_biasBound_config <- function(X, 
                                   Y = NULL,
                                   h = 0.09, 
                                   alpha = 0.05, 
                                   resol = 100,
                                   xi_lb = NULL, 
                                   xi_ub = NULL, 
                                   methods_get_xi = "Schennach",
                                   kernel.fun = "Schennach2004", 
                                   if_approx_kernel = TRUE, 
                                   kernel.resol = 1000) {
  
  # Create kernel functions
  kernel_functions <- create_kernel_functions(
    kernel.fun = kernel.fun,
    if_approx_kernel = if_approx_kernel,
    kernel.resol = kernel.resol,
    X = X,
    h = h
  )
  
  # Get xi interval if not specified
  if (is.null(xi_lb) || is.null(xi_ub)) {
    xi_interval <- get_xi_interval(X = X, methods = methods_get_xi)
    xi_lb <- xi_interval$xi_lb
    xi_ub <- xi_interval$xi_ub
  } else {
    xi_interval <- list(xi_lb = xi_lb, xi_ub = xi_ub)
  }
  
  # Get A and r estimates
  est_Ar <- get_est_Ar(X = X, xi_interval = xi_interval)
  
  # Get b1x estimation
  b1x <- get_est_b1x(
    X = X, 
    h = h, 
    est_Ar = est_Ar, 
    inf_k_ft = kernel_functions$kernel_ft
  )
  
  # Get B and byx if Y is provided
  if (!is.null(Y)) {
    est_B <- get_est_B(Y = Y)
    byx <- get_est_byx(
      Y = Y, 
      X = X, 
      h = h, 
      est_Ar = est_Ar, 
      est_B = est_B, 
      inf_k_ft = kernel_functions$kernel_ft
    )
  } else {
    est_B <- NULL
    byx <- NULL
  }
  
  # Create and return config object
  config <- list(
    X = X,
    Y = Y,
    h = h,
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
    est_Ar = est_Ar,
    b1x = b1x,
    est_B = est_B,
    byx = byx
  )
  
  return(config)
} 