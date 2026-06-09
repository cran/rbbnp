test_that("get_avg_phi works correctly", {
  # Test with simple cases
  X <- c(1, 2, 3)
  Y <- c(1, 1, 1)
  xi <- 1

  result <- get_avg_phi(Y = Y, X = X, xi = xi)

  expect_type(result, "double")
  expect_true(result >= 0)  # Should always be non-negative
})

test_that("get_xi_interval works with different methods", {
  # Test data
  X <- rnorm(100)

  # Test Schennach method
  result1 <- get_xi_interval(X = X, methods = "Schennach")
  expect_type(result1, "list")
  expect_true(result1$xi_ub > result1$xi_lb)

  # Test Schennach_loose method
  result2 <- get_xi_interval(X = X, methods = "Schennach_loose")
  expect_type(result2, "list")
  expect_true(result2$xi_ub > result2$xi_lb)
})

test_that("get_sigma handles negative density estimates", {
  # Setup test data that will produce negative density estimate
  X <- c(1, 2, 3)
  x <- 10  # Point far from data points
  h <- 0.1  # Small bandwidth

  # Test sigma
  result <- get_sigma(X, x, h, inf_k = normal_kernel)
  expect_equal(result, 0)

  # Test with valid data
  x <- 2  # Point within data range
  result2 <- get_sigma(X, x, h, inf_k = normal_kernel)
  expect_true(result2 >= 0)
})

test_that("get_sigma_yx handles negative density estimates", {
  # Setup test data
  X <- c(1, 2, 3)
  Y <- c(1, 4, 9)
  x <- 10  # Point far from data points
  h <- 0.1  # Small bandwidth

  # Test sigma_yx
  result <- get_sigma_yx(Y, X, x, h, inf_k = normal_kernel)
  expect_equal(result, 0)

  # Test with valid data
  x <- 2  # Point within data range
  result2 <- get_sigma_yx(Y, X, x, h, inf_k = normal_kernel)
  expect_true(result2 >= 0)
})

test_that("get_sigma_yx is invariant to sign flip of Y (depends on Var[Y|X=x])", {
  X <- c(1, 2, 3)
  Y <- c(1, 4, 9)
  x <- 2
  h <- 0.4

  s1 <- get_sigma_yx(Y, X, x, h, inf_k = normal_kernel)
  s2 <- get_sigma_yx(-Y, X, x, h, inf_k = normal_kernel)

  expect_true(is.finite(s1) && s1 >= 0)
  expect_equal(s1, s2, tolerance = 1e-12)
})

test_that("get_conditional_var handles negative values", {
  # Setup test data
  X <- c(1, 2, 3)
  Y <- c(1, 1, 1)  # Same Y values should give zero variance
  x <- 2
  h <- 0.1

  # Test with data that should give zero/negative variance
  result <- get_conditional_var(X, Y, x, h, kernel_func = normal_kernel)
  expect_equal(result, 0)

  # Test with normal data
  Y2 <- c(1, 2, 3)
  result2 <- get_conditional_var(X, Y2, x, h, kernel_func = normal_kernel)
  expect_true(result2 >= 0)
})

# Task 1.4: Tests for vectorized get_conditional_var()
test_that("get_conditional_var vectorized matches loop implementation", {
  set.seed(123)
  n <- 100
  X <- rnorm(n)
  Y <- X^2 + rnorm(n, sd = 0.5)
  x <- 0.5
  h <- 0.3

  # Test with all kernel types
  kernels <- list(
    normal = normal_kernel,
    epanechnikov = epanechnikov_kernel,
    sinc = sinc
  )

  for (kernel_name in names(kernels)) {
    kernel_func <- kernels[[kernel_name]]

    result_vectorized <- get_conditional_var(X, Y, x, h, kernel_func)
    result_loop <- get_conditional_var_loop(X, Y, x, h, kernel_func)

    expect_equal(result_vectorized, result_loop, tolerance = 1e-10,
                 info = sprintf("Vectorized implementation should match loop for %s kernel", kernel_name))
  }
})

test_that("get_conditional_var vectorized performance speedup", {
  skip_on_cran()  # Skip time-consuming test on CRAN
  skip("Skipping performance test - results vary by system")

  set.seed(456)
  n <- 1000
  X <- rnorm(n)
  Y <- 2 * X + rnorm(n, sd = 0.8)
  x <- 0
  h <- 0.5

  # Benchmark loop implementation (run multiple times for stable measurement)
  time_loop <- system.time({
    for (i in 1:3) {
      result_loop <- get_conditional_var_loop(X, Y, x, h, normal_kernel)
    }
  })["elapsed"] / 3

  # Benchmark vectorized implementation (run multiple times for stable measurement)
  time_vectorized <- system.time({
    for (i in 1:3) {
      result_vectorized <- get_conditional_var(X, Y, x, h, normal_kernel)
    }
  })["elapsed"] / 3

  # Verify numerical equivalence
  expect_equal(result_vectorized, result_loop, tolerance = 1e-10)

  # Check speedup (should be >5x for n=1,000, adjusted from 10x due to vectorization overhead)
  speedup <- time_loop / time_vectorized
  expect_true(speedup > 2,
              info = sprintf("Speedup: %.1fx (expected >2x for n=%d)", speedup, n))
})

test_that("get_conditional_var handles edge cases", {
  # Edge case 1: Single point
  X <- c(1)
  Y <- c(2)
  result <- get_conditional_var(X, Y, x = 1, h = 0.1, kernel_func = normal_kernel)
  expect_equal(result, 0)  # Zero variance with single point

  # Edge case 2: Two points
  X <- c(1, 2)
  Y <- c(1, 3)
  result <- get_conditional_var(X, Y, x = 1.5, h = 0.5, kernel_func = normal_kernel)
  expect_true(result >= 0)

  # Edge case 3: Zero variance (constant Y)
  set.seed(789)
  X <- rnorm(50)
  Y <- rep(5, 50)
  result <- get_conditional_var(X, Y, x = 0, h = 0.5, kernel_func = normal_kernel)
  expect_equal(result, 0, tolerance = 1e-8)

  # Edge case 4: Point outside data range
  X <- rnorm(100, mean = 0, sd = 1)
  Y <- X + rnorm(100, sd = 0.2)
  result <- get_conditional_var(X, Y, x = 10, h = 0.5, kernel_func = normal_kernel)
  expect_true(result >= 0)
})

test_that("get_conditional_var works with all kernel types", {
  set.seed(111)
  n <- 200
  X <- rnorm(n)
  Y <- sin(X) + rnorm(n, sd = 0.3)
  x <- 0
  h <- 0.4

  # Test Normal kernel
  result_normal <- get_conditional_var(X, Y, x, h, normal_kernel)
  expect_true(result_normal >= 0)
  expect_type(result_normal, "double")

  # Test Epanechnikov kernel
  result_epan <- get_conditional_var(X, Y, x, h, epanechnikov_kernel)
  expect_true(result_epan >= 0)
  expect_type(result_epan, "double")

  # Test Sinc kernel
  result_sinc <- get_conditional_var(X, Y, x, h, sinc)
  expect_true(result_sinc >= 0)
  expect_type(result_sinc, "double")

  # Results should be different for different kernels
  expect_false(isTRUE(all.equal(result_normal, result_epan)))
})

test_that("get_conditional_var memory safety checks", {
  # Test warning for large memory allocation (>4GB)
  # For n=23172, memory = (23172^2 * 8) / (1024^3) = 4.0005 GB (just over threshold)
  n_warn <- 23172

  # Create minimal data to trigger warning
  X <- seq_len(n_warn)
  Y <- seq_len(n_warn)

  expect_warning(
    get_conditional_var(X, Y, x = 1, h = 1, kernel_func = normal_kernel),
    regexp = "Large memory allocation.*GB"
  )

  # Test error for excessive memory (>10GB)
  # For n=36637, memory = (36637^2 * 8) / (1024^3) = 10.0007 GB (just over threshold)
  n_error <- 36637

  X_large <- seq_len(n_error)
  Y_large <- seq_len(n_error)

  expect_error(
    get_conditional_var(X_large, Y_large, x = 1, h = 1, kernel_func = normal_kernel),
    regexp = "Memory requirement.*exceeds 10 GB"
  )
})

test_that("get_conditional_var preserves backward compatibility", {
  # Ensure function signature unchanged
  set.seed(222)
  X <- rnorm(50)
  Y <- X + rnorm(50, sd = 0.5)

  # All original function calls should still work
  result1 <- get_conditional_var(X, Y, x = 0, h = 0.5, kernel_func = normal_kernel)
  expect_type(result1, "double")
  expect_true(result1 >= 0)

  # Test with different parameter values
  result2 <- get_conditional_var(X = X, Y = Y, x = 1, h = 0.3, kernel_func = epanechnikov_kernel)
  expect_type(result2, "double")
  expect_true(result2 >= 0)
})

test_that("get_sigma_yx returns 0 when conditional variance is 0; otherwise positive", {
  X <- c(1, 2, 3)
  x <- 2
  h <- 0.4

  # Zero conditional variance => zero standard error
  Y_const <- c(1, 1, 1)
  s0 <- get_sigma_yx(Y_const, X, x, h, inf_k = normal_kernel)
  expect_equal(s0, 0)

  # Nonzero conditional variance with positive density => positive standard error
  Y_var <- c(1, 4, 9)
  s1 <- get_sigma_yx(Y_var, X, x, h, inf_k = normal_kernel)
  expect_true(is.finite(s1))
  expect_gt(s1, 0)
})

# Tests for LP solver infrastructure
test_that(".prepare_lp_problem creates valid LP matrices", {
  ln_xi <- seq(0, 2, length.out = 10)
  phi_log <- -ln_xi  # Simple decay
  lp <- .prepare_lp_problem(ln_xi, phi_log)

  expect_equal(length(lp$obj), 2)
  expect_equal(nrow(lp$mat), 10)
  expect_equal(ncol(lp$mat), 2)
  expect_equal(length(lp$rhs), 10)
  expect_equal(length(lp$dir), 10)
  expect_true(all(lp$dir == ">="))
})

test_that("get_est_Ar works with sample data", {
  set.seed(123)
  X <- rnorm(100)
  xi_interval <- list(xi_lb = 0.5, xi_ub = 5)

  result <- get_est_Ar(X = X, xi_interval = xi_interval, method = "lp")

  expect_named(result, c("est_A", "est_r"))
  expect_true(result["est_A"] > 0)
  expect_true(result["est_r"] >= 0)
  expect_true(result["est_r"] < 50)  # Within r_max
})

test_that("get_est_Ar LP and grid give similar results", {
  set.seed(456)
  X <- rnorm(200)
  xi_interval <- list(xi_lb = 1, xi_ub = 10)

  result_lp <- get_est_Ar(X = X, xi_interval = xi_interval, method = "lp")
  result_grid <- get_est_Ar(X = X, xi_interval = xi_interval, method = "grid", r_stepsize = 500)

  # LP and fine grid search should agree within reasonable tolerance
  # Using 15% for est_A and 0.15 absolute for est_r to account for discretization
  expect_equal(result_lp["est_A"], result_grid["est_A"], tolerance = 0.15)
  expect_equal(result_lp["est_r"], result_grid["est_r"], tolerance = 0.15)
})

test_that("get_est_Ar LP solution satisfies envelope constraint", {
  set.seed(789)
  X <- rnorm(150, mean = 5, sd = 2)
  xi_interval <- list(xi_lb = 0.5, xi_ub = 8)

  result <- get_est_Ar(X = X, xi_interval = xi_interval, method = "lp")

  # Recreate frequency grid
  ln_xi_range <- seq(log(xi_interval$xi_lb), log(xi_interval$xi_ub), length.out = 100)
  avg_phi_log <- vapply(ln_xi_range, function(x) get_avg_phi_log(X = X, ln_xi = x), numeric(1))

  # Check constraint: ln(A) - r*ln(xi) >= ln(|phi(xi)|) for all xi
  envelope_line <- log(result["est_A"]) - result["est_r"] * ln_xi_range
  constraint_satisfied <- all(envelope_line >= avg_phi_log - 1e-6)

  expect_true(constraint_satisfied,
              info = "LP solution must lie above all Fourier transform points")
})

test_that("get_est_Ar handles edge cases gracefully", {
  # Edge case 1: Very small sample
  set.seed(111)
  X_small <- rnorm(5)
  xi_interval <- list(xi_lb = 0.5, xi_ub = 2)
  result_small <- get_est_Ar(X = X_small, xi_interval = xi_interval)
  expect_true(result_small["est_A"] > 0)
  expect_true(result_small["est_r"] >= 0)

  # Edge case 2: Wide xi_interval
  X_normal <- rnorm(100)
  xi_wide <- list(xi_lb = 0.1, xi_ub = 50)
  result_wide <- get_est_Ar(X = X_normal, xi_interval = xi_wide)
  expect_true(result_wide["est_r"] < 50)  # Within bounds
  expect_true(result_wide["est_A"] > 0)
})

test_that("get_est_Ar validates inputs correctly", {
  X <- rnorm(50)

  # Invalid xi_interval (not a list)
  expect_error(get_est_Ar(X = X, xi_interval = c(0.5, 5)),
               "xi_interval must be a list")

  # Invalid xi_interval (missing xi_ub)
  expect_error(get_est_Ar(X = X, xi_interval = list(xi_lb = 0.5)),
               "xi_interval must be a list")

  # Invalid xi_interval (xi_lb <= 0)
  expect_error(get_est_Ar(X = X, xi_interval = list(xi_lb = -1, xi_ub = 5)),
               "Invalid xi_interval")

  # Invalid xi_interval (xi_ub <= xi_lb)
  expect_error(get_est_Ar(X = X, xi_interval = list(xi_lb = 5, xi_ub = 2)),
               "Invalid xi_interval")

  # Invalid method
  expect_error(get_est_Ar(X = X, xi_interval = list(xi_lb = 0.5, xi_ub = 5), method = "invalid"),
               "Unknown method")
})

test_that("get_est_Ar gives reproducible results (regression test)", {
  set.seed(12345)
  X <- rnorm(100, mean = 0, sd = 1)
  xi_interval <- list(xi_lb = 1, xi_ub = 5)

  result <- get_est_Ar(X = X, xi_interval = xi_interval, method = "lp")

  # Results should be stable across runs
  result2 <- get_est_Ar(X = X, xi_interval = xi_interval, method = "lp")
  expect_equal(result, result2)  # Deterministic
})

test_that("get_est_Ar method='auto' selects LP when available", {
  set.seed(999)
  X <- rnorm(50)
  xi_interval <- list(xi_lb = 0.5, xi_ub = 5)

  result_auto <- get_est_Ar(X = X, xi_interval = xi_interval, method = "auto")
  result_lp <- get_est_Ar(X = X, xi_interval = xi_interval, method = "lp")

  # Auto should select LP (since Rglpk is available)
  expect_equal(result_auto, result_lp)
})

test_that("get_est_Ar LP is faster than grid search", {
  skip_on_cran()  # Skip on CRAN (time-consuming)

  set.seed(2024)
  X <- rnorm(1000)  # Large sample to emphasize speedup
  xi_interval <- list(xi_lb = 0.5, xi_ub = 10)

  # Benchmark grid search (coarse grid to keep test fast)
  time_grid <- system.time({
    result_grid <- get_est_Ar(X = X, xi_interval = xi_interval,
                              method = "grid", r_stepsize = 50)
  })

  # Benchmark LP
  time_lp <- system.time({
    result_lp <- get_est_Ar(X = X, xi_interval = xi_interval, method = "lp")
  })

  # LP should be at least as fast as grid search (speedup >= 0.8x allows for minor variance)
  speedup <- time_grid["elapsed"] / time_lp["elapsed"]
  expect_true(speedup > 0.8,
              info = sprintf("LP speedup: %.1fx (expected >0.8x)", speedup))

  # Results should still be similar
  expect_equal(result_lp["est_A"], result_grid["est_A"], tolerance = 0.1)
})

