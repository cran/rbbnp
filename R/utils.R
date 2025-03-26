#' Compute Sample Average of Fourier Transform Magnitude
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param xi A single numerical value representing the frequency at which the Fourier transform
#'    is computed.
#'
#' @return Returns the sample estimation of expected Fourier transform at frequency `xi`.
#'
get_avg_phi <- function(Y = 1, X, xi) {
  return(
    (mean(Y * cos(xi * X))^2 + mean(Y * sin(xi * X))^2)^0.5
  )
}

#' Compute log sample average of fourier transform and get mod
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param ln_xi A single numerical value representing the log frequency at which the Fourier transform
#'    is computed.
#'
#' @return Returns the log sample estimation of expected Fourier transform at frequency `xi`.
get_avg_phi_log <- function(Y = 1, X, ln_xi) {
  return(
    log(get_avg_phi(Y = Y, X = X, xi = exp(ln_xi)))
  )
}

#' get the estimation of Vy
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#'
get_est_vy <- function(Y) {
  mean(Y^2)^0.5
}

#' get xi interval
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param methods A character string indicating the method to use for calculating the
#'        xi interval. Supported methods are "Schennach" and "Schennach_loose".
#'        Defaults to "Schennach".
#'
#' @return A list containing the lower (`xi_lb`) and upper (`xi_ub`) bounds of the xi interval.
#'
#' @details
#' The "Schennach" method computes the xi interval by performing a test based on the
#' Schennach's theorem, adjusting the upper bound `xi_ub` if the test condition is met.
#' The "Schennach_loose" method provides a looser calculation of the xi interval without
#' performing the Schennach's test.
#'
get_xi_interval <- function(Y = 1, X, methods = "Schennach") {
  if (methods == "Schennach") {
    n <- length(X)
    avg_dphi <- 7^0.5 * get_est_vy(Y) * n^(-0.5) * log(n)^0.5
    min_avg_phi <- avg_dphi * log(n) # the minimum |\hat_Φ| can pass the test

    # lower bound of xi
    avg_v_x <- var(X)^0.5
    xi_lb <- 1 / avg_v_x

    xi_ub <- n^0.25 / avg_v_x

    # Define the default xi_n by theorem 2
    xi_range <- seq(xi_lb, xi_ub, length.out = log(xi_ub / xi_lb) * 500)

    # calculate the Fourier transform vector
    avg_phi <- purrr::map_dbl(xi_range, function(u) get_avg_phi(Y = Y, X = X, xi = u))

    # refine the xi_n by finding the max xi in [xi_lb, xi_ub]
    # which let avg_phi larger than min_avg_phi
    xi_n_feasible <- xi_range[which(avg_phi > min_avg_phi)]

    if (length(xi_n_feasible) == 0) {
      warning("No feasible xi_n can pass the Schennach's test in [xi_lb, xi_ub]")
      xi_n <- NA
    } else {
      xi_n <- max(xi_n_feasible)
    }

    # if Schennach test generate non empty results, then assign it to xi_ub
    # otherwise directly use the xi_ub
    if (!is.na(xi_n)) {
      xi_ub <- xi_n
    }

    return(list(xi_lb = xi_lb, xi_ub = xi_ub))
  }

  if (methods == "Schennach_loose") {
    n <- length(X)
    avg_dphi <- 7^0.5 * get_est_vy(Y) * n^(-0.5) * log(n)^0.5
    min_avg_phi <- avg_dphi * log(n) # the minimum |\hat_Φ| can pass the test

    # lower bound of xi
    avg_v_x <- var(X)^0.5
    xi_lb <- 1 / avg_v_x

    xi_ub <- n^0.25 / avg_v_x
    return(list(xi_lb = xi_lb, xi_ub = xi_ub))
  }
}

#' get the estimation of A and r
#'
#' This function estimates the parameters A and r by optimizing an objective function
#' over a specified range of frequency values and r values.
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param xi_interval A list with elements `xi_lb` and `xi_ub` representing the lower
#'        and upper bounds of the frequency interval.
#' @param r_stepsize An integer value representing the number of steps in the r range.
#'        This controls the granularity of the estimation. Higher values lead to finer
#'        granularity but increase computation time.
#'
#' @return A named vector with elements `est_A` and `est_r` representing the estimated
#'         values of A and r, respectively.
#'
#' @details
#' The function internally defines a range for the natural logarithm of frequency values (`ln_xi_range`)
#' and a range for the parameter `r` (`r_range`). It then defines an optimization function `optim_ln_A`
#' to minimize the integral of a given function over the `ln_xi_range`. The actual estimation is done by
#' finding the `r` and `A` value that minimizes the the area of the line \eqn{\ln A - r \ln \xi} under the constraint that the line should not go below the Fourier transform curve.
#'
get_est_Ar <- function(Y = 1, X, xi_interval, r_stepsize = 150) {
  # get the interval of xi
  ln_xi_lb <- log(xi_interval$xi_lb)
  ln_xi_ub <- log(xi_interval$xi_ub)
  ln_xi_range <- seq(ln_xi_lb, ln_xi_ub, length.out = (ln_xi_ub - ln_xi_lb) * 200)

  # Precompute get_avg_phi_log values for all ln_xi
  avg_phi_log_values <- vapply(ln_xi_range, function(x) get_avg_phi_log(X = X, ln_xi = x), numeric(1))

  # define the range of r
  r_range <- tan(seq(-pi / 2, pi / 2, length.out = r_stepsize))
  r_range <- r_range[-c(1, length(r_range))]

  # Initial ln_A value (using the first precomputed value)
  ln_A_init <- avg_phi_log_values[1]

  # Define the Optimize lnA for a given ln_xi with precomputed values
  optim_ln_A <- function(r) {
    # Define the objective function to minimize (integral)
    objective_function <- function(ln_A) {
      ln_xi_lb <- ln_xi_range[1]
      ln_xi_n <- ln_xi_range[length(ln_xi_range)]
      # the integration can be viewed as the area of a rectangular
      (ln_xi_n - ln_xi_lb) * (ln_A - r * (ln_xi_n + ln_xi_lb) / 2) # length*height
    }

    # Compute constraint values using precomputed phi_log values
    ln_phi_c <- ln_A_init - r * ln_xi_range - avg_phi_log_values

    diff_ln_phi <- min(ln_phi_c)

    # get the optimal ln_A
    ln_A <- ln_A_init - diff_ln_phi

    obj <- objective_function(ln_A)

    return(c(ln_A = ln_A, obj = obj))
  }

  # Create a results matrix
  res <- matrix(0, nrow = length(r_range), ncol = 2)
  colnames(res) <- c("ln_A", "obj")

  # Evaluate for each r value
  for (i in seq_along(r_range)) {
    res[i,] <- optim_ln_A(r_range[i])
  }

  # Convert to data frame
  res <- as.data.frame(res)

  # find the optimal value of A and r
  est_id <- which.min(res$obj)
  est_A <- exp(res$ln_A[est_id])
  est_r <- r_range[est_id]

  return(c(est_A = est_A, est_r = est_r))
}

#' get the estimation of B
#' @param Y A numerical vector representing the sample data of variable Y.
#'
#' @return The mean of the absolute values of the elements in Y, representing the estimated value of \eqn{B}.
#'
get_est_B <- function(Y) {
  est_B <- mean(abs(Y))
  return(est_B)
}

#' Estimation of bias b1x
#'
#' Computes the bias estimate for given parameters.
#'
#' @param X A numerical vector representing the sample data of variable X.
#' @param ... Additional arguments passed to other methods.
#' @return A scalar representing the bias b1x estimate.
get_est_b1x <- function(X, ...) {

  # Integration expression of b1x
  b1x_int <- function(xi, h, est_Ar = get_est_Ar(X = X), inf_k_ft = W_kernel_ft) {
    Mod(1 - inf_k_ft(xi * h)) * apply(cbind(1, est_Ar[1] * abs(xi)^(-est_Ar[2])), FUN = min, MARGIN = 1)
  }

  pracma::quadgk(b1x_int, a = -99999, b = 100000, ...) / (2 * pi)
}


#' Estimation of bias byx
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param ... Additional arguments passed to other methods.
#' @return A scalar representing the bias byx estimate.
get_est_byx <- function(Y, X, ...) {

  # Integration expression of byx
  byx_int <- function(xi, h, est_Ar = get_est_Ar(Y = Y, X = X), est_B = get_est_B(Y = Y),
                      inf_k_ft = W_kernel_ft) {
    Mod(1 - inf_k_ft(xi * h)) * apply(cbind(est_B, est_Ar[1] * abs(xi)^(-est_Ar[2])), FUN = min, MARGIN = 1)
  }

  pracma::quadgk(byx_int, a = -99999, b = 100000, ...) / (2 * pi)
}

#' Kernel point estimation
#'
#' Computes the point estimate using the specified kernel function.
#'
#' @param X A numerical vector of sample data.
#' @param x A scalar representing the point where the density is estimated.
#' @param h A scalar bandwidth parameter.
#' @param inf_k Kernel function used for the computation.
#' @return A scalar representing the kernel density estimate at point x.
get_avg_f1x <- function(X, x, h, inf_k) {
  mean(inf_k(u = (x - X) / h) / h)
}

#' Kernel point estimation
#'
#' Computes the point estimate using the specified kernel function.
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param x A scalar representing the point where the density is estimated.
#' @param h A scalar bandwidth parameter.
#' @param inf_k Kernel function used for the computation.
#' @return A scalar representing the kernel density estimate at point x.
get_avg_fyx <- function(Y, X, x, h, inf_k) {
  mean(Y * inf_k(u = (x - X) / h) / h)
}

#' Kernel Regression function
#'
#' @param X A numerical vector representing the sample data of variable X.
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param x The point at which the regression function is to be estimated.
#' @param h A bandwidth parameter that determines the weight assigned to each observation in X.
#' @param kernel_func A function that computes the weight of each observation based on its distance to x.
#'
#' @return Returns a scalar representing the estimated value of the regression function at the point x.
#'
kernel_reg <- function(X, Y, x, h, kernel_func) {
  weights <- kernel_func((x - X) / h)
  return(sum(weights * Y) / sum(weights))
}

#' get the conditional variance of Y on X for given x
#'
#' @param X A numerical vector representing the sample data of variable X.
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param x The specific point at which the conditional variance is to be calculated.
#' @param h A bandwidth parameter used in the kernel function for smoothing.
#' @param kernel_func A kernel function used to weigh observations in the neighborhood of point x.
#'
#' @return Returns a non-negative scalar representing the estimated conditional variance of Y given X at the point x.
#'         Returns 0 if the computed variance is negative.
#'
get_conditional_var <- function(X, Y, x, h, kernel_func) {
  # get the conditional mean of Y on X
  conditional_mean <- sapply(X, function(x) {
    kernel_reg(X, Y, x, h, kernel_func)
  })

  # Compute residuals
  residuals <- Y - conditional_mean

  # Compute conditional variance
  conditional_var <- kernel_reg(X, residuals^2, x = x, h, kernel_func = kernel_func)
  
  # Return 0 if variance is negative, otherwise return the computed variance
  return(max(0, conditional_var))
}

#' Estimation of sigma
#'
#' Computes the sigma estimate for given parameters.
#'
#' @param X A numerical vector of sample data.
#' @param x A scalar representing the point where the density is estimated.
#' @param h A scalar bandwidth parameter.
#' @param inf_k Kernel function used for the computation.
#' @return A scalar representing the sigma estimate at point x. Returns 0 if the density estimate is negative.
get_sigma <- function(X, x, h, inf_k) {
  n <- length(X)
  f1x <- get_avg_f1x(X, x, h, inf_k = inf_k)
  
  # Return 0 if density estimate is negative or 0
  if (f1x <= 0) {
    return(0)
  }
  
  z <- f1x / (n * h) * pracma::integral(function(x) {
    inf_k(x)^2
  }, -999, +1000)
  
  return(sqrt(z))
}

#' Estimation of sigma_yx
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param x The specific point at which sigma_yx is to be estimated.
#' @param h A bandwidth parameter used in the kernel function for smoothing.
#' @param inf_k A kernel function used to weigh observations in the neighborhood of point x.
#'
#' @return Returns a scalar representing the estimated value of sigma_yx at the point x. 
#'         Returns 0 if either fyx or conditional variance is negative.
get_sigma_yx <- function(Y, X, x, h, inf_k) {
  n <- length(X)
  var_yx <- get_conditional_var(Y = Y, X = X, x = x, h = h, kernel_func = inf_k)
  fyx <- get_avg_fyx(Y = Y, X = X, x = x, h = h, inf_k = inf_k)
  
  # Return 0 if either fyx is negative or var_yx is 0 (indicating it was negative)
  if (fyx <= 0 || var_yx <= 0) {
    return(0)
  }
  
  z <- var_yx * fyx / (n * h) * pracma::integral(function(x) {
    inf_k(x)^2
  }, -999, +1000)
  
  return(sqrt(z))
}

#' Approximation Function for Intensive Calculations
#'
#' This function provides a lookup-based approximation for calculations that are computationally intensive.
#' Once computed, it stores the results in an environment and uses linear interpolation for new data points
#' to speed up subsequent computations.
#'
#' @param u A vector of values where the function should be evaluated.
#' @param u_lb Lower bound for the precomputed range. Defaults to -10.
#' @param u_ub Upper bound for the precomputed range. Defaults to 10.
#' @param resol The resolution or number of sample points in the precomputed range. Defaults to 1000.
#' @param fun A function for which the approximation is computed. Defaults to the `W` function.
#'
#' @return A vector of approximated function values corresponding to `u`.
#'
#' @details
#' The `fun_approx` function works by initially creating a lookup table of function values based on
#' the range specified by `u_lb` and `u_ub` and the resolution `resol`. This precomputation only happens once
#' for a given set of parameters (`u_lb`, `u_ub`, `resol`, and `fun`). Subsequent calls to `fun_approx` with the
#' same parameters use the lookup table to find the closest precomputed points to the requested `u` values
#' and then return an interpolated result.
#'
#' Linear interpolation is used between the two closest precomputed points in the lookup table. This
#' ensures a smooth approximation for values in between sample points.
#'
#' This function is especially useful for computationally intensive functions where recalculating
#' function values is expensive or time-consuming. By using a combination of precomputation and
#' interpolation, `fun_approx` provides a balance between accuracy and speed.
#'
fun_approx <- (function() {
  # Initialize an environment to store the interpolated function, bounds, and resolution
  storage_env <- new.env()

  # The inner approximation function
  function(u, u_lb = -100, u_ub = 100, resol = 1000, fun = W_kernel) {
    # Check if it's the first call or if parameters have changed
    if (!exists("local_interpolated_fun", envir = storage_env) ||
      storage_env$local_lb != u_lb ||
      storage_env$local_ub != u_ub ||
      storage_env$local_resol != resol) {
      sample_points <- seq(u_lb, u_ub, length.out = resol)
      precomputed_W <- fun(sample_points)
      interpolated_fun <- approxfun(sample_points, precomputed_W, method = "linear")

      # Store values in the environment
      storage_env$local_interpolated_fun <- interpolated_fun
      storage_env$local_lb <- u_lb
      storage_env$local_ub <- u_ub
      storage_env$local_resol <- resol
    }

    # Use the stored interpolated function
    res <- storage_env$local_interpolated_fun(u)

    # if the value exceed the range, NA will generated
    # replace those NA as 0 since most of the kernel function is close to 0 at tail
    return(ifelse(is.na(res), 0, res))
  }
})()

