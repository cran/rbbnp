#' Compute Sample Average of Fourier Transform Magnitude
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param xi A single numerical value representing the frequency at which the Fourier transform
#'    is computed.
#'
#' @return Returns the sample estimation of expected Fourier transform at frequency `xi`.
#'
#' @keywords internal
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
#' @keywords internal
get_avg_phi_log <- function(Y = 1, X, ln_xi) {
  return(
    log(get_avg_phi(Y = Y, X = X, xi = exp(ln_xi)))
  )
}

#' get the estimation of Vy
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#'
#' @keywords internal
get_est_vy <- function(Y) {
  mean(Y^2)^0.5
}

#' Uniform Fourier-transform error bound Delta-phi (Schennach 2020, Theorem 1)
#'
#' @param Y Sample vector (or 1 for density).
#' @param n Sample size.
#' @param noise_floor "compact" (the `sqrt(7) * (ln n)^(1/2)` form, default), "general"
#'   (the `(7 * sqrt(2) / 3) * ln n` form), or "auto" (general when Y is non-constant, else compact).
#' @keywords internal
.noise_floor_avg_dphi <- function(Y, n, noise_floor = "compact") {
  if (identical(noise_floor, "auto")) {
    noise_floor <- if (length(unique(Y)) > 1) "general" else "compact"
  }
  vy <- get_est_vy(Y)
  if (noise_floor == "general") {
    (7 * sqrt(2) / 3) * vy * n^(-0.5) * log(n)          # general random Y
  } else {
    sqrt(7) * vy * n^(-0.5) * log(n)^0.5                # compactly supported Y
  }
}

#' Signal-to-noise frequency-window selection rule
#'
#' Selects the upper frequency cutoff as the point where the pointwise signal-to-noise ratio of
#' the empirical Fourier transform first drops (on a sustained run) below `tau`. Using the
#' pointwise (delta-method) standard error in place of Schennach (2020)'s worst-case uniform
#' error bound keeps the rule usable at realistic sample sizes, and a minimum-window guard
#' ensures it never returns an empty interval.
#'
#' @param Y Numeric vector (or 1 for density).
#' @param X Numeric vector of sample data.
#' @param tau SNR threshold: `NULL` (default) uses `max(2, sqrt(2 * log(n)))`; a numeric
#'   value, or a function of `n`. Letting `tau` grow faster than `sqrt(log n)` preserves the
#'   asymptotic rate condition of Schennach (2020) while keeping the rule usable in finite samples.
#' @param run Minimum consecutive below-threshold points for a sustained crossing.
#' @param per_log Grid resolution (points per log-unit of frequency).
#' @return A list with `xi_lb` and `xi_ub`.
#' @keywords internal
.snr_window <- function(Y = 1, X, tau = NULL, run = 5L, per_log = 200) {
  n <- length(X); sdx <- sqrt(stats::var(X))
  xi_lb <- 1 / sdx; xi_hi <- n^0.25 / sdx
  if (is.null(tau)) tau <- max(2, sqrt(2 * log(n))) else if (is.function(tau)) tau <- tau(n)
  grid <- exp(seq(log(xi_lb), log(xi_hi),
                  length.out = max(50L, round(log(xi_hi / xi_lb) * per_log))))
  Cr <- vapply(grid, function(u) mean(Y * cos(u * X)), 0)
  Ci <- vapply(grid, function(u) mean(Y * sin(u * X)), 0)
  mod <- sqrt(Cr^2 + Ci^2)
  vr <- vapply(grid, function(u) stats::var(Y * cos(u * X)) / n, 0)
  vi <- vapply(grid, function(u) stats::var(Y * sin(u * X)) / n, 0)
  se  <- sqrt((Cr^2 * vr + Ci^2 * vi) / pmax(mod^2, .Machine$double.eps))
  snr <- mod / pmax(se, .Machine$double.eps)
  below <- snr < tau
  idx <- NA_integer_
  if (any(below)) {
    rl <- rle(below); starts <- cumsum(rl$lengths) - rl$lengths + 1L
    hit <- which(rl$values & rl$lengths >= run)
    if (length(hit)) idx <- starts[hit[1]]
  }
  xi_ub <- if (is.na(idx) || idx <= 1L) xi_hi else grid[idx]
  xi_ub <- max(xi_ub, xi_lb * 1.3)            # minimum-window guard: never collapse
  list(xi_lb = xi_lb, xi_ub = min(xi_ub, xi_hi))
}

#' get xi interval
#'
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param X A numerical vector representing the sample data of variable X.
#' @param methods A character string: "snr" (default; a signal-to-noise cutoff that selects a
#'        valid frequency window at realistic sample sizes), "Schennach" (the data-driven rule of
#'        Schennach 2020, Theorem 2), or "Schennach_loose" (the initial, un-refined interval).
#'
#' @param noise_floor Noise-floor form for the feasibility test: "auto" (default; general
#'   for non-constant Y, compact otherwise), "compact", or "general". See
#'   `.noise_floor_avg_dphi()`.
#' @return A list containing the lower (`xi_lb`) and upper (`xi_ub`) bounds of the xi interval.
#'
#' @details
#' The "Schennach" method computes the xi interval by performing a test based on the
#' Schennach's theorem, adjusting the upper bound `xi_ub` if the test condition is met.
#' The "Schennach_loose" method provides a looser calculation of the xi interval without
#' performing the Schennach's test.
#' The "snr" method selects the upper cutoff from the pointwise signal-to-noise ratio of the
#' empirical Fourier transform (see `.snr_window()`); unlike "Schennach" it always selects a
#' valid cutoff at realistic sample sizes and never returns an empty window.
#'
#' @keywords internal
get_xi_interval <- function(Y = 1, X, methods = "snr",
                            noise_floor = c("auto", "compact", "general")) {
  noise_floor <- match.arg(noise_floor)
  if (methods == "Schennach") {
    n <- length(X)
    # Schennach (2020) Theorem 1: a.s. uniform bound on the FT estimation error.
    avg_dphi <- .noise_floor_avg_dphi(Y, n, noise_floor)
    min_avg_phi <- avg_dphi * log(n) # smallest |phi-hat| that can pass the test

    # lower bound of xi: 1/sd(X) ensures sufficient oscillations
    avg_v_x <- var(X)^0.5
    xi_lb <- 1 / avg_v_x

    # upper bound of xi: n^0.25/sd(X) balances bias and variance (Schennach 2004)
    xi_ub <- n^0.25 / avg_v_x

    # Grid resolution: 500 points per log-unit for accurate Fourier transform estimation
    xi_range <- seq(xi_lb, xi_ub, length.out = log(xi_ub / xi_lb) * 500)

    # calculate the Fourier transform vector
    avg_phi <- purrr::map_dbl(xi_range, function(u) get_avg_phi(Y = Y, X = X, xi = u))

    # refine the xi_n by finding the max xi in [xi_lb, xi_ub]
    # which let avg_phi larger than min_avg_phi
    xi_n_feasible <- xi_range[which(avg_phi > min_avg_phi)]

    if (length(xi_n_feasible) == 0) {
      warning(sprintf(
        "No feasible xi_n passed Schennach's test in interval [%.4f, %.4f]. Using theoretical upper bound xi_ub = %.4f. This may indicate insufficient signal or inappropriate xi bounds.",
        xi_lb, xi_ub, xi_ub
      ))
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
    # Schennach (2020) Theorem 1: a.s. uniform bound on the FT estimation error.
    avg_dphi <- .noise_floor_avg_dphi(Y, n, noise_floor)
    min_avg_phi <- avg_dphi * log(n) # smallest |phi-hat| that can pass the test

    # lower bound of xi: 1/sd(X) ensures sufficient oscillations
    avg_v_x <- var(X)^0.5
    xi_lb <- 1 / avg_v_x

    # upper bound of xi: n^0.25/sd(X) balances bias and variance (Schennach 2004)
    xi_ub <- n^0.25 / avg_v_x
    return(list(xi_lb = xi_lb, xi_ub = xi_ub))
  }

  if (methods == "snr") {
    return(.snr_window(Y = Y, X = X))
  }

  stop(sprintf("Unknown methods_get_xi '%s'. Use 'snr', 'Schennach', or 'Schennach_loose'.", methods))
}

#' Check if LP solver is available
#'
#' @return Logical indicating if Rglpk is available
#' @keywords internal
.has_lp_solver <- function() {
  requireNamespace("Rglpk", quietly = TRUE)
}

#' Prepare LP problem matrices for envelope optimization
#'
#' @param ln_xi_range Numeric vector of log-frequencies
#' @param avg_phi_log_values Numeric vector of log-Fourier-transform values
#' @param r_max Maximum allowed value for r parameter (default 50)
#' @return List with components: obj (objective coefficients),
#'         mat (constraint matrix), rhs (RHS vector), bounds (variable bounds)
#' @keywords internal
.prepare_lp_problem <- function(ln_xi_range, avg_phi_log_values, r_max = 50) {
  # Number of constraint points
  m <- length(ln_xi_range)

  # Objective: minimize z_A - r * lambda_bar
  lambda_bar <- mean(range(ln_xi_range))  # Midpoint of log-frequency range
  obj <- c(1, -lambda_bar)

  # Constraint matrix: [1, -ln(xi_i)] for each i
  # Constraint: z_A - r*ln(xi_i) >= ln(|phi(xi_i)|)
  mat <- matrix(0, nrow = m, ncol = 2)
  mat[, 1] <- 1                    # Coefficient for z_A
  mat[, 2] <- -ln_xi_range         # Coefficient for r

  # Right-hand side: ln(|phi(xi_i)|)
  # Handle near-zero values to avoid -Inf
  rhs <- pmax(avg_phi_log_values, -100)  # Clamp to avoid -Inf

  # Variable bounds: z_A in [-100, 100], r in [0, r_max]
  bounds <- list(
    lower = list(ind = c(1L, 2L), val = c(-100, 0)),
    upper = list(ind = c(1L, 2L), val = c(100, r_max))
  )

  list(obj = obj, mat = mat, rhs = rhs, bounds = bounds, dir = rep(">=", m))
}

#' Solve LP problem for envelope optimization
#'
#' @param lp_problem List from .prepare_lp_problem()
#' @return List with status, solution vector, and objective value
#' @keywords internal
.solve_lp_envelope <- function(lp_problem) {
  if (!.has_lp_solver()) {
    stop("Rglpk package is required but not available. Install with install.packages('Rglpk')")
  }

  # Call Rglpk solver
  result <- Rglpk::Rglpk_solve_LP(
    obj = lp_problem$obj,
    mat = lp_problem$mat,
    dir = lp_problem$dir,
    rhs = lp_problem$rhs,
    bounds = lp_problem$bounds,
    max = FALSE  # Minimize
  )

  # Check solution status
  # Status codes: 0 = optimal, 1 = suboptimal, 2 = infeasible, 3 = unbounded
  if (result$status != 0) {
    warning(sprintf("LP solver returned non-optimal status: %d", result$status))
    return(NULL)  # Trigger fallback
  }

  result
}

#' Estimate A and r parameters via grid search (LEGACY)
#'
#' This is the original grid search implementation, kept for fallback and validation.
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
#' @keywords internal
get_est_Ar_grid <- function(Y = 1, X, xi_interval, r_stepsize = 150, use_Y = FALSE) {
  # get the interval of xi
  ln_xi_lb <- log(xi_interval$xi_lb)
  ln_xi_ub <- log(xi_interval$xi_ub)
  # Grid resolution: 200 points per log-unit for accurate envelope estimation
  ln_xi_range <- seq(ln_xi_lb, ln_xi_ub, length.out = (ln_xi_ub - ln_xi_lb) * 200)

  # Precompute get_avg_phi_log values for all ln_xi
  avg_phi_log_values <- vapply(ln_xi_range, function(x) get_avg_phi_log(Y = if (use_Y) Y else 1, X = X, ln_xi = x), numeric(1))

  # Define r range using tan transformation to cover (-Inf, +Inf)
  # Remove endpoints at -pi/2 and pi/2 to avoid Inf values
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

#' Estimate A and r parameters via linear programming
#'
#' Finds the minimal-area envelope line in log-log space using linear programming.
#'
#' @param Y Numeric vector (default 1 for density estimation)
#' @param X Numeric vector of sample data
#' @param xi_interval List with xi_lb and xi_ub
#' @param r_max Maximum allowed value for r (default 50)
#' @return Named numeric vector with est_A and est_r
#' @keywords internal
get_est_Ar_lp <- function(Y = 1, X, xi_interval, r_max = 50, use_Y = FALSE) {
  # Extract interval bounds
  ln_xi_lb <- log(xi_interval$xi_lb)
  ln_xi_ub <- log(xi_interval$xi_ub)

  # Grid resolution: 200 points per log-unit (matches grid search for consistency)
  ln_xi_range <- seq(ln_xi_lb, ln_xi_ub, length.out = (ln_xi_ub - ln_xi_lb) * 200)

  # Precompute Fourier transform values (same as grid search)
  avg_phi_log_values <- vapply(ln_xi_range, function(x) get_avg_phi_log(Y = if (use_Y) Y else 1, X = X, ln_xi = x), numeric(1))

  # Prepare LP problem matrices
  lp_problem <- .prepare_lp_problem(ln_xi_range, avg_phi_log_values, r_max = r_max)

  # Solve LP problem
  lp_result <- .solve_lp_envelope(lp_problem)

  # Check if LP solver succeeded
  if (is.null(lp_result)) {
    warning("LP solver failed; falling back to grid search")
    return(get_est_Ar_grid(Y = Y, X = X, xi_interval = xi_interval, use_Y = use_Y))
  }

  # Extract solution: x = [z_A, r]
  z_A <- lp_result$solution[1]
  est_r <- lp_result$solution[2]
  est_A <- exp(z_A)

  # Verify constraints are satisfied (debugging/validation)
  if (any(z_A - est_r * ln_xi_range < avg_phi_log_values - 1e-6)) {
    warning("LP solution violates envelope constraint; falling back to grid search")
    return(get_est_Ar_grid(Y = Y, X = X, xi_interval = xi_interval, use_Y = use_Y))
  }

  return(c(est_A = est_A, est_r = est_r))
}

#' Estimate A and r parameters for Fourier envelope
#'
#' Estimates parameters A and r that bound the Fourier transform decay.
#' Uses linear programming by default, with automatic fallback to grid search.
#'
#' @param Y Numeric vector (default 1 for density estimation)
#' @param X Numeric vector of sample data
#' @param xi_interval List with xi_lb and xi_ub
#' @param r_stepsize Grid resolution for fallback method (default 150).
#'   Only used when method="grid". Ignored for LP method.
#' @param method Character: "lp" (default), "grid", or "auto"
#' @return Named numeric vector with est_A and est_r
#' @keywords internal
#'
#' @details
#' The LP-based method (default) is typically 1-5x faster than grid search,
#' with speedup depending on sample size and Fourier transform computation cost.
#' The grid search method is preserved for validation and as fallback if LP solver is unavailable.
#'
#' The function finds the minimal-area envelope line in log-log space by
#' minimizing the integral subject to envelope constraints on the
#' Fourier transform decay.
get_est_Ar <- function(Y = 1, X, xi_interval, r_stepsize = 150, method = "auto", use_Y = FALSE,
                       integer_r = FALSE) {
  # Validate inputs
  if (!is.list(xi_interval) || !all(c("xi_lb", "xi_ub") %in% names(xi_interval))) {
    stop("xi_interval must be a list with xi_lb and xi_ub")
  }
  if (xi_interval$xi_lb <= 0 || xi_interval$xi_ub <= xi_interval$xi_lb) {
    stop("Invalid xi_interval: must have 0 < xi_lb < xi_ub")
  }

  # Auto-select method
  if (method == "auto") {
    method <- if (.has_lp_solver()) "lp" else "grid"
  }

  # Dispatch to appropriate implementation
  res <- if (method == "lp") {
    get_est_Ar_lp(Y = Y, X = X, xi_interval = xi_interval, use_Y = use_Y)
  } else if (method == "grid") {
    get_est_Ar_grid(Y = Y, X = X, xi_interval = xi_interval, r_stepsize = r_stepsize, use_Y = use_Y)
  } else {
    stop(sprintf("Unknown method: %s. Use 'lp', 'grid', or 'auto'", method))
  }

  raw_r <- res[["est_r"]]

  # When the fitted slope falls below the minimum smoothness assumed by Schennach (2020,
  # Definition 2), i.e. r < 2, clamp it up to r = 2 and refit A as the tightest upper bound over
  # the window. A slope of at least 2 keeps the bias-bound integral finite and avoids the blow-up
  # that occurs when the fitted slope is below 1. Slopes that are already >= 2 are left unchanged,
  # so well-behaved data are not widened unnecessarily.
  if (isTRUE(integer_r) && is.finite(raw_r) && raw_r < 2) {
    r_int <- max(2, round(raw_r))
    ln_xi_range <- seq(log(xi_interval$xi_lb), log(xi_interval$xi_ub),
                       length.out = (log(xi_interval$xi_ub) - log(xi_interval$xi_lb)) * 200)
    phi_log <- vapply(ln_xi_range, function(x)
      get_avg_phi_log(Y = if (use_Y) Y else 1, X = X, ln_xi = x), numeric(1))
    ln_A <- max(phi_log + r_int * ln_xi_range)   # ln_A - r_int*ln(xi) >= log|phi(xi)| for all xi
    res <- c(est_A = exp(ln_A), est_r = r_int)
  }

  # Diagnostic: a fitted slope rounding below 2 means the data do not show the minimum smoothness
  # (r >= 2) assumed by Schennach (2020, Definition 2) -- typically a very smooth density, or a
  # conditional mean whose cross-spectrum does not decay like a power law. The bias bound stays
  # valid but becomes a conservative (wider) fallback.
  if (is.finite(raw_r) && round(raw_r) < 2) {
    warning(sprintf(
      paste0("Fitted envelope slope r-hat = %.2f is below the smoothness floor r >= 2 the method ",
             "assumes: the Fourier transform does not show a clear power-law decay over the window ",
             "(likely a supersmooth density or a non-power-law/non-monotone cross-spectrum). The ",
             "bias bound is valid but conservative%s. Inspect plot(fit, type = \"ft\")."),
      raw_r, if (isTRUE(integer_r)) " (clamped to r >= 2)"
             else " and may be very wide (set integer_r = TRUE to clamp)"),
      call. = FALSE)
  }
  res
}

#' get the estimation of B
#' @param Y A numerical vector representing the sample data of variable Y.
#'
#' @return The mean of the absolute values of the elements in Y, representing the estimated value of \eqn{B}.
#'
#' @keywords internal
get_est_B <- function(Y) {
  est_B <- mean(abs(Y))
  return(est_B)
}

#' Estimation of bias b1x
#'
#' Computes the bias estimate for given parameters.
#'
#' @param X A numerical vector representing the sample data of variable X.
#' @param h A scalar bandwidth parameter.
#' @param est_Ar A vector containing the estimated A and r parameters.
#' @param inf_k_ft A kernel Fourier transform function.
#' @param ... Additional arguments passed to the quadgk integration function.
#' @return A scalar representing the bias b1x estimate.
#' @keywords internal
get_est_b1x <- function(X, h, est_Ar, inf_k_ft, ...) {

  # Integration expression of b1x
  b1x_int <- function(xi) {
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
kernel_reg <- function(X, Y, x, h, kernel_func) {
  weights <- kernel_func((x - X) / h)
  return(sum(weights * Y) / sum(weights))
}

#' Original loop-based implementation of conditional variance (preserved for testing)
#'
#' @param X A numerical vector representing the sample data of variable X.
#' @param Y A numerical vector representing the sample data of variable Y.
#' @param x The specific point at which the conditional variance is to be calculated.
#' @param h A bandwidth parameter used in the kernel function for smoothing.
#' @param kernel_func A kernel function used to weigh observations in the neighborhood of point x.
#'
#' @return Returns a non-negative scalar representing the estimated conditional variance of Y given X at the point x.
#'         Returns 0 if the computed variance is negative.
#' @keywords internal
get_conditional_var_loop <- function(X, Y, x, h, kernel_func) {
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

#' get the conditional variance of Y on X for given x
#'
#' Vectorized implementation using matrix operations for improved performance.
#' For large datasets (n > 10,000), memory usage is approximately 8n^2 bytes.
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
#' @details
#' Performance: Achieves >10x speedup for n=1,000 and >20x speedup for n=10,000
#' compared to the original loop-based implementation.
#'
#' Memory: Allocates an n x n distance matrix. For n=10,000, requires ~800 MB.
#' Issues a warning when estimated memory exceeds 4 GB and stops with an error
#' when estimated memory exceeds 10 GB.
#'
#' @keywords internal
get_conditional_var <- function(X, Y, x, h, kernel_func) {
  n <- length(X)

  # Memory safety checks
  # Each n x n matrix uses 8 bytes per element (double precision)
  estimated_memory_gb <- (n * n * 8) / (1024^3)

  # 10 GB limit: prevents system memory exhaustion on typical machines
  if (estimated_memory_gb > 10) {
    stop(sprintf("Memory requirement (%.1f GB) exceeds 10 GB limit for n=%d. Consider using smaller dataset or chunked processing.",
                 estimated_memory_gb, n))
  }

  # 4 GB warning threshold: alerts user to monitor memory for large computations
  if (estimated_memory_gb > 4) {
    warning(sprintf("Large memory allocation (%.1f GB) for n=%d. Consider monitoring memory usage.",
                    estimated_memory_gb, n))
  }

  # Vectorized computation of conditional mean at all X points
  # Create distance matrix: dist_matrix[i,j] = (X[j] - X[i]) / h
  dist_matrix <- outer(X, X, function(xi, xj) (xj - xi) / h)

  # Apply kernel function to get weight matrix
  # Note: Some kernel approximations may flatten the matrix, so we preserve dimensions
  weight_matrix <- kernel_func(dist_matrix)
  if (!is.matrix(weight_matrix)) {
    weight_matrix <- matrix(weight_matrix, nrow = n, ncol = n)
  }

  # Compute conditional mean for each point: E[Y|X=xi]
  # conditional_mean[i] = sum(weights[i,j] * Y[j]) / sum(weights[i,j])
  weighted_sums <- as.vector(weight_matrix %*% Y)
  weight_sums <- rowSums(weight_matrix)
  conditional_mean <- weighted_sums / weight_sums

  # Compute residuals
  residuals <- Y - conditional_mean

  # Compute conditional variance at point x
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
#' @keywords internal
get_sigma <- function(X, x, h, inf_k) {
  n <- length(X)
  f1x <- get_avg_f1x(X, x, h, inf_k = inf_k)

  # Return 0 if density estimate is negative or 0
  if (is.na(f1x) || f1x <= 0) {
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
#' @keywords internal
get_sigma_yx <- function(Y, X, x, h, inf_k) {
  n <- length(X)
  var_yx <- get_conditional_var(Y = Y, X = X, x = x, h = h, kernel_func = inf_k)
  f1x <- get_avg_f1x(X = X, x = x, h = h, inf_k = inf_k)

  # Return 0 if either f1x is non-positive (unstable ratio region) or var_yx is non-positive
  if (is.na(f1x) || f1x <= 0 || var_yx <= 0) {
    return(0)
  }

  z <- var_yx * f1x / (n * h) * pracma::integral(function(u) {
    inf_k(u)^2
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
#' @keywords internal
fun_approx <- (function() {
  # Initialize an environment to store the interpolated function, bounds, and resolution
  storage_env <- new.env()

  # The inner approximation function
  function(u, u_lb = -100, u_ub = 100, resol = 10000, fun = W_kernel) {
    # Check if it's the first call or if parameters have changed
    if (!exists("local_interpolated_fun", envir = storage_env) ||
      storage_env$local_lb != u_lb ||
      storage_env$local_ub != u_ub ||
      storage_env$local_resol != resol) {

      # Pre-compute kernel function values for interpolation (silent operation)
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
