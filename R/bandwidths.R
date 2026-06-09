#' @importFrom stats median sd IQR fft
NULL

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

# Internal helper: Linear binning of data onto regular grid
# Returns list with counts vector and grid parameters for FFT-based density estimation
.bin_data_linear <- function(X, grid_size = 512, extend_factor = 3, h_ref = NULL) {
  n <- length(X)

  # Round grid_size up to next power of 2 for FFT efficiency
  M <- 2^ceiling(log2(grid_size))

  # Compute grid bounds with extension
  # Use h_ref * extend_factor beyond data range to avoid boundary effects
  if (is.null(h_ref)) {
    # If no h_ref provided, use a simple estimate based on data range
    h_ref <- diff(range(X)) / 10
  }

  x_min <- min(X) - extend_factor * h_ref
  x_max <- max(X) + extend_factor * h_ref

  # Create regular grid
  grid <- seq(x_min, x_max, length.out = M)
  delta <- (x_max - x_min) / (M - 1)

  # Initialize counts vector
  counts <- rep(0, M)

  # Bin each data point with linear interpolation weights
  for (i in seq_along(X)) {
    # Find position in grid
    pos <- (X[i] - x_min) / delta + 1  # +1 for 1-based indexing

    # Handle edge cases: points outside grid go to boundary bins
    if (pos < 1) {
      counts[1] <- counts[1] + 1
    } else if (pos > M) {
      counts[M] <- counts[M] + 1
    } else {
      # Linear binning: distribute mass between two nearest bins
      j <- floor(pos)
      w_right <- pos - j  # Weight for right bin
      w_left <- 1 - w_right  # Weight for left bin

      # Accumulate counts
      if (j >= 1 && j <= M) {
        counts[j] <- counts[j] + w_left
      }
      if (j + 1 >= 1 && j + 1 <= M) {
        counts[j + 1] <- counts[j + 1] + w_right
      }
    }
  }

  # Validate: sum(counts) should equal n (within numerical tolerance)
  total_count <- sum(counts)
  if (abs(total_count - n) > 1e-6) {
    warning(sprintf("Binning error: sum(counts) = %.6f, expected n = %d", total_count, n))
  }

  return(list(
    counts = counts,
    grid = grid,
    delta = delta,
    M = M
  ))
}

# Internal helper: Compute kernel density via FFT convolution
# X_binned: output from .bin_data_linear()
# h: bandwidth
# kernel_func: kernel function
# Returns: density values on grid
.fft_convolve_density <- function(X_binned, h, kernel_func) {
  n <- sum(X_binned$counts)  # Total number of observations
  M <- X_binned$M
  grid <- X_binned$grid

  # Evaluate kernel on grid scaled by h
  K_grid <- kernel_func(grid / h) / h

  # Apply FFT to kernel and binned counts
  K_fft <- stats::fft(K_grid)
  counts_fft <- stats::fft(X_binned$counts / n)

  # Multiply in frequency domain (convolution theorem)
  density_fft <- K_fft * counts_fft

  # Inverse FFT to get density
  density <- stats::fft(density_fft, inverse = TRUE) / M

  # Extract real part (imaginary part should be negligible)
  density_real <- Re(density)

  # Validate: check imaginary part is negligible
  max_imag <- max(abs(Im(density)))
  if (max_imag > 1e-8) {
    warning(sprintf("Large imaginary component in FFT result: max(|Im|) = %.2e", max_imag))
  }

  # Clamp near-zero values to exactly zero (avoid negative densities from rounding)
  density_real[density_real < 0 & density_real > -1e-10] <- 0

  return(density_real)
}

# Internal helper: Compute LSCV score for given bandwidth via FFT
# X_binned: binned data
# h: bandwidth
# kernel_func: kernel function
# X_original: original data points (for leave-one-out correction)
# Returns: LSCV score (scalar)
.lscv_score_fft <- function(X_binned, h, kernel_func, X_original) {
  n <- length(X_original)

  # Compute density via FFT convolution
  density_grid <- .fft_convolve_density(X_binned, h, kernel_func)

  # term1: integral of squared density = sum(f^2) * delta
  term1 <- sum(density_grid^2) * X_binned$delta

  # term2: leave-one-out correction
  # Map original X to grid indices for interpolation
  x_min <- X_binned$grid[1]
  delta <- X_binned$delta

  # Find density at original data points using linear interpolation
  pos <- (X_original - x_min) / delta + 1  # Position in grid (1-based)

  # Linear interpolation for each data point
  f_at_X <- numeric(n)
  for (i in seq_along(X_original)) {
    j <- floor(pos[i])

    # Handle edge cases
    if (j < 1) {
      f_at_X[i] <- density_grid[1]
    } else if (j >= X_binned$M) {
      f_at_X[i] <- density_grid[X_binned$M]
    } else {
      # Linear interpolation between grid points
      w <- pos[i] - j
      f_at_X[i] <- (1 - w) * density_grid[j] + w * density_grid[j + 1]
    }
  }

  # Leave-one-out correction: subtract self-contribution K(0)/(nh)
  self_contrib <- kernel_func(0) / (n * h)
  term2 <- (2 / n) * sum(f_at_X - self_contrib)

  # LSCV score
  cv_score <- term1 - term2

  return(cv_score)
}

#' Cross-Validation for Bandwidth Selection
#'
#' Implements least-squares cross-validation for bandwidth selection with any kernel function.
#' Uses FFT-based algorithm for n >= 100 (fast, O(m log m) complexity) and exact computation
#' for n < 100 (accurate). The FFT method bins data onto a regular grid and computes the
#' LSCV objective via convolution in the frequency domain.
#'
#' @param X A numerical vector of sample data.
#' @param h_grid A numerical vector of bandwidth values to evaluate. If NULL (default), a grid is
#'        automatically generated based on the range and distribution of the data.
#' @param kernel_func The kernel function to use for cross-validation.
#' @param kernel_type A string identifying the kernel type, used only for reference bandwidth.
#' @param grid_size Number of grid points for FFT-based evaluation (used when n >= 100).
#'        Default is 512. Larger values increase accuracy but reduce speed. Automatically
#'        rounded up to the next power of 2.
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
    h_grid <- seq(h_ref * 0.2, h_ref * 2, length.out = 100)
  }

  # Dispatch: use exact method for small n, FFT method for large n
  if (n < 100) {
    # Use exact implementation (preserve for small n accuracy)
    return(.cv_bandwidth_exact(X, h_grid, kernel_func))
  } else {
    # Use FFT-based fast implementation
    return(.cv_bandwidth_fft(X, h_grid, kernel_func, grid_size))
  }
}

# Internal function: Exact cross-validation (direct O(n^2) implementation)
# Used for n < 100 where accuracy is prioritized over speed
.cv_bandwidth_exact <- function(X, h_grid, kernel_func) {
  n <- length(X)

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

# Internal function: FFT-based cross-validation (O(m log m) implementation)
# Used for n >= 100 where speed is prioritized
.cv_bandwidth_fft <- function(X, h_grid, kernel_func, grid_size = 512) {
  n <- length(X)

  # Estimate reference bandwidth for grid extension
  h_ref <- median(h_grid)  # Use median of h_grid as typical bandwidth

  # Bin data onto regular grid (only once for all bandwidths)
  X_binned <- .bin_data_linear(X, grid_size = grid_size, extend_factor = 5, h_ref = h_ref)

  # Compute LSCV score for each bandwidth via FFT
  cv_scores <- vapply(h_grid, function(h) {
    .lscv_score_fft(X_binned, h, kernel_func, X)
  }, numeric(1))

  # Return bandwidth that minimizes CV score
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




