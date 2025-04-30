#' Silverman's Rule of Thumb for Bandwidth Selection
#'
#' Implements Silverman's rule of thumb for selecting an optimal bandwidth in kernel density estimation.
#'
#' @param X A numerical vector of sample data.
#' @param kernel_type A string identifying the kernel type.
#'
#' @return A scalar representing the optimal bandwidth.
#'
#' @export
#'
#' @examples
#' # Generate sample data
#' X <- rnorm(100)
#' # Get optimal bandwidth using Silverman's rule
#' h_opt <- silverman_bandwidth(X, kernel_type = "normal")
silverman_bandwidth <- function(X, kernel_type = "normal") {
  n <- length(X)
  sd_X <- sd(X)
  iqr_X <- IQR(X)

  # Use min of SD and normalized IQR to make robust against outliers
  spread <- min(sd_X, iqr_X/1.34)

  # Standard Silverman factor for Gaussian kernel
  factor <- 0.9

  # Adjust factor based on kernel type
  if (kernel_type == "epanechnikov") {
    factor <- 2.34 * factor  # Adjustment for Epanechnikov kernel
  } else if (kernel_type %in% c("Schennach2004", "sinc")) {
    # For infinite-order kernels, the adjustment is less standardized
    factor <- 1.2 * factor
  }

  # Calculate optimal bandwidth using Silverman's rule
  h_opt <- factor * spread * n^(-1/5)

  return(h_opt)
}

#' Cross-Validation for Bandwidth Selection
#'
#' Implements least-squares cross-validation for bandwidth selection with any kernel function.
#' Uses the self-convolution approach for accurate estimation of the integral term.
#'
#' @param X A numerical vector of sample data.
#' @param h_grid A numerical vector of bandwidth values to evaluate. If NULL (default), a grid is
#'        automatically generated based on the range and distribution of the data.
#' @param kernel_func The kernel function to use for cross-validation.
#' @param kernel_type A string identifying the kernel type, used only for reference bandwidth.
#' @param grid_size Number of grid points for evaluation. Default is 512.
#'
#' @return A scalar representing the optimal bandwidth that minimizes the cross-validation score.
#'
#' @export
#'
#' @examples
#' # Generate sample data
#' X <- rnorm(100)
#' # Get optimal bandwidth using cross-validation with a normal kernel
#' kernel_functions <- create_kernel_functions("normal")
#' h_opt <- cv_bandwidth(X, kernel_func = kernel_functions$kernel,
#'                      kernel_type = kernel_functions$kernel_type)
cv_bandwidth <- function(X, h_grid = NULL, kernel_func, kernel_type = "normal", grid_size = 512) {
  n <- length(X)

  # Check if X has sufficient elements for cross-validation
  if (n < 3) {
    stop("At least 3 data points are needed for cross-validation")
  }

  # Silverman's rule of thumb for bandwidth reference
  if (is.null(h_grid)) {
    h_ref <- silverman_bandwidth(X, kernel_type)
    message(paste0("h_ref = ", round(h_ref, 6)))
    h_grid <- seq(h_ref * 0.2, h_ref * 2, length.out = 100)
  }

  precompute_self_convolution <- function(kernel, h, resol = 1000) {
    u_grid <- seq(-10*h, 10*h, length.out = resol)  # assume the kernel is supported in [-10h, 10h]
    k_conv <- sapply(u_grid, function(u) {
      pracma::integral(function(v) kernel(v) * kernel(u - v),
                xmin = -10*h, xmax = 10*h)
    })
    approxfun(u_grid, k_conv, rule = 2)
  }
  self_conv <- precompute_self_convolution(kernel = kernel_func, h = 5)

  # Pre-compute matrices for efficiency
  X_mat <- outer(X, X, "-")

  # CV for each h
  cv_scores <- sapply(h_grid, function(h) {
    # Calculate term1 using self-convolution method
    # Compute pairwise differences scaled by bandwidth
    scaled_diffs <- X_mat / h

    # Apply self-convolution to each pairwise difference
    convolution_values <- self_conv(as.vector(scaled_diffs))

    # Compute term1 as mean of all self-convolutions divided by bandwidth
    term1 <- mean(convolution_values) / h

    # Second term: average of leave-one-out density estimates
    # Apply kernel function to normalized differences and preserve matrix structure
    K_mat <- matrix(kernel_func(X_mat / h), nrow = n, ncol = n)

    # Set diagonal elements to 0 for leave-one-out
    diag(K_mat) <- 0

    # Second term: average of leave-one-out density estimates
    loo_sums <- rowSums(K_mat)
    term2 <- 2 * mean(loo_sums / ((n - 1) * h))

    return(term1 - term2)
  })

  # Return the bandwidth that minimizes the CV score
  h_opt <- h_grid[which.min(cv_scores)]
  return(h_opt)
}

#' Select Optimal Bandwidth
#'
#' Selects an optimal bandwidth using the specified method.
#'
#' @param X A numerical vector of sample data.
#' @param Y Optional. A numerical vector of sample data for conditional expectation estimation.
#' @param method A string specifying the bandwidth selection method. Options are "cv" for
#'        cross-validation and "silverman" for Silverman's rule of thumb. Defaults to "cv".
#' @param kernel.fun A string specifying the kernel type. Options include "normal", "epanechnikov",
#'        "Schennach2004", and "sinc".
#' @param if_approx_kernel Logical. If TRUE, uses approximations for the kernel function.
#' @param kernel.resol The resolution for kernel function approximation.
#'
#' @return A scalar representing the optimal bandwidth.
#'
#' @export
#'
#' @examples
#' # Generate sample data
#' X <- rnorm(100)
#' # Get optimal bandwidth using cross-validation with normal kernel
#' h_opt <- select_bandwidth(X, method = "cv", kernel.fun = "normal")
#' # Get optimal bandwidth using Silverman's rule with Schennach kernel
#' h_opt <- select_bandwidth(X, method = "silverman", kernel.fun = "Schennach2004")
select_bandwidth <- function(X, Y = NULL, method = "cv", kernel.fun = "normal",
                           if_approx_kernel = TRUE, kernel.resol = 1000) {

  # Ensure kernel.fun is a string
  if (!is.character(kernel.fun)) {
    stop("kernel.fun must be a string specifying the kernel type (e.g., 'normal', 'Schennach2004')")
  }

  # Create kernel functions from the string
  message(paste0("create kernel functions for ", kernel.fun))
  kernel_functions <- create_kernel_functions(
    kernel.fun = kernel.fun,
    if_approx_kernel = if_approx_kernel,
    kernel.resol = kernel.resol
  )

  # Extract the kernel function and type from the list
  kernel_func <- kernel_functions$kernel
  kernel_type <- kernel_functions$kernel_type

  # For conditional expectation with Y data, use the same bandwidth as for X
  # This is a simplification - ideally we'd use a specific method for conditional expectation
  if (method == "cv") {
    return(cv_bandwidth(X, kernel_func = kernel_func, kernel_type = kernel_type))
  } else if (method == "silverman") {
    return(silverman_bandwidth(X, kernel_type = kernel_type))
  } else {
    stop("Unsupported bandwidth selection method. Use 'cv' or 'silverman'.")
  }
}




