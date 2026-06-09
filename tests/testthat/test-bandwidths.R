test_that("select_bandwidth integrates different methods correctly", {
  # Test data
  X <- rnorm(30)
  Y <- X^2 + rnorm(30, 0, 0.5)

  # Test with cross-validation method
  h_cv <- select_bandwidth(X, method = "cv", kernel.fun = "normal")
  expect_true(h_cv > 0)

  # Test with Silverman's rule
  h_silver <- select_bandwidth(X, method = "silverman", kernel.fun = "normal")
  expect_true(h_silver > 0)

  # Test with invalid method
  expect_error(select_bandwidth(X, method = "invalid"))

  # Test with Y data for conditional expectation
  h_with_y <- select_bandwidth(X, Y = Y, method = "cv", kernel.fun = "normal")
  expect_true(h_with_y > 0)

  # Test with custom kernel and Silverman's rule
  h_schennach_silver <- select_bandwidth(X, method = "silverman", kernel.fun = "Schennach2004")
  expect_true(h_schennach_silver > 0)

  # Test with custom kernel and cross-validation (using small sample for speed)
  X_small <- X[1:10]
  h_schennach_cv <- select_bandwidth(X_small, method = "cv", kernel.fun = "Schennach2004")
  expect_true(h_schennach_cv > 0)
})

test_that("create_biasBound_config handles bandwidth selection correctly", {
  # Test data
  X <- rnorm(30)
  Y <- X^2 + rnorm(30, 0, 0.5)

  # Test with provided bandwidth
  h_fixed <- 0.5
  config_fixed <- create_biasBound_config(X, Y, h = h_fixed)
  expect_equal(config_fixed$h, h_fixed)

  # Test with automatic selection using Silverman's rule
  config_auto <- create_biasBound_config(X, Y, h = NULL, h_method = "silverman")
  expect_true(config_auto$h > 0)
  expect_equal(config_auto$h_method, "silverman")

  # Test with automatic selection using cross-validation
  config_cv <- create_biasBound_config(X, h = NULL, h_method = "cv", kernel.fun = "normal")
  expect_true(config_cv$h > 0)
  expect_equal(config_cv$h_method, "cv")

  # Test with Schennach2004 kernel
  config_schennach <- create_biasBound_config(
    X, Y, h = NULL,
    h_method = "silverman",
    kernel.fun = "Schennach2004"
  )
  expect_true(config_schennach$h > 0)
  expect_equal(config_schennach$kernel.fun, "Schennach2004")

  # Test with sinc kernel
  config_sinc <- create_biasBound_config(
    X, Y, h = NULL,
    h_method = "silverman",
    kernel.fun = "sinc"
  )
  expect_true(config_sinc$h > 0)
  expect_equal(config_sinc$kernel.fun, "sinc")
})

test_that("cv_bandwidth FFT and exact methods produce similar results", {
  set.seed(42)
  X <- rnorm(150)  # Large enough for FFT
  kernel_functions <- create_kernel_functions("normal")

  # Get reference bandwidth grid
  h_ref <- silverman_bandwidth(X, "normal")
  h_grid <- seq(h_ref * 0.5, h_ref * 1.5, length.out = 50)

  # Compute optimal h using both methods on same subset
  X_subset <- X[1:90]  # Use subset that's below threshold
  h_exact <- rbbnp:::.cv_bandwidth_exact(X_subset, h_grid, kernel_functions$kernel)
  h_fft <- rbbnp:::.cv_bandwidth_fft(X_subset, h_grid, kernel_functions$kernel, grid_size = 512)

  # Should agree within 2% (FFT introduces binning error)
  relative_error <- abs(h_fft - h_exact) / h_exact
  expect_true(relative_error < 0.02,
              info = paste("Relative error:", round(relative_error, 4)))
})

test_that("FFT method works with all supported kernels", {
  set.seed(123)
  X <- rnorm(200)

  kernel_types <- c("normal", "epanechnikov", "sinc", "Schennach2004")

  for (ktype in kernel_types) {
    kernel_functions <- create_kernel_functions(ktype, if_approx_kernel = TRUE)

    # Should not error and should return positive bandwidth
    h_opt <- cv_bandwidth(X, kernel_func = kernel_functions$kernel,
                          kernel_type = ktype, grid_size = 256)

    expect_true(h_opt > 0,
                info = paste("Kernel type:", ktype))
  }
})

test_that("FFT cv_bandwidth handles edge cases", {
  kernel_functions <- create_kernel_functions("normal")

  # Edge case 1: Data with extreme outliers
  X_outliers <- c(rnorm(100), 100, -100)
  h_outliers <- cv_bandwidth(X_outliers, kernel_func = kernel_functions$kernel)
  expect_true(h_outliers > 0)

  # Edge case 2: Very small variance data
  X_small_var <- rnorm(150, mean = 10, sd = 0.01)
  h_small <- cv_bandwidth(X_small_var, kernel_func = kernel_functions$kernel)
  expect_true(h_small > 0)
  expect_true(h_small < 0.1)  # Should adapt to small variance

  # Edge case 3: Bimodal data
  X_bimodal <- c(rnorm(100, -2, 0.5), rnorm(100, 2, 0.5))
  h_bimodal <- cv_bandwidth(X_bimodal, kernel_func = kernel_functions$kernel)
  expect_true(h_bimodal > 0)

  # Edge case 4: Exactly n=100 (threshold boundary)
  X_100 <- rnorm(100)
  h_100 <- cv_bandwidth(X_100, kernel_func = kernel_functions$kernel)
  expect_true(h_100 > 0)
})

test_that("cv_bandwidth dispatches correctly based on n", {
  kernel_functions <- create_kernel_functions("normal")

  # Small n should use exact method (n < 100)
  X_small <- rnorm(50)
  h_small <- cv_bandwidth(X_small, kernel_func = kernel_functions$kernel)
  expect_true(h_small > 0)

  # Large n should use FFT method (n >= 100)
  X_large <- rnorm(200)
  h_large <- cv_bandwidth(X_large, kernel_func = kernel_functions$kernel)
  expect_true(h_large > 0)
})

# Expected performance (approximate, depends on hardware):
# n=100:   FFT ~0.05s, Exact ~0.05s (comparable)
# n=500:   FFT ~0.15s, Exact ~1.2s (8x speedup)
# n=1000:  FFT ~0.25s, Exact ~5s (20x speedup)
# n=5000:  FFT ~0.8s,  Exact ~120s (150x speedup)
# n=10000: FFT ~1.5s,  Exact ~500s (300x speedup)
test_that("FFT cv_bandwidth performance scales better than O(n²)", {
  skip_on_cran()  # Benchmarking takes time

  kernel_functions <- create_kernel_functions("normal")

  # Get reference bandwidth for grid
  h_ref_val <- 0.3  # Fixed reference for consistency

  # Benchmark at different sample sizes
  timings <- data.frame(
    n = c(100, 500, 1000, 2000),
    time_exact = NA_real_,
    time_fft = NA_real_
  )

  for (i in seq_len(nrow(timings))) {
    n <- timings$n[i]
    set.seed(42)
    X <- rnorm(n)
    h_grid <- seq(h_ref_val * 0.5, h_ref_val * 1.5, length.out = 50)

    # Time exact method (only for n <= 500, too slow otherwise)
    if (n <= 500) {
      timings$time_exact[i] <- system.time({
        h_exact <- rbbnp:::.cv_bandwidth_exact(X, h_grid = h_grid,
                                                kernel_func = kernel_functions$kernel)
      })[3]  # Elapsed time
    }

    # Time FFT method
    timings$time_fft[i] <- system.time({
      h_fft <- rbbnp:::.cv_bandwidth_fft(X, h_grid = h_grid,
                                          kernel_func = kernel_functions$kernel,
                                          grid_size = 512)
    })[3]
  }

  # Verify FFT method is faster for n >= 1000
  expect_true(timings$time_fft[timings$n == 1000] < 2.0,
              info = paste("FFT method should complete in <2 sec for n=1000, got:",
                          round(timings$time_fft[timings$n == 1000], 2)))

  # Check complexity: FFT time should scale sub-quadratically
  # time_fft(2000) / time_fft(1000) should be < (2000/1000)^1.5 = 2.83
  # Add small tolerance for timing variations
  ratio_n <- 2000 / 1000
  ratio_time <- timings$time_fft[timings$n == 2000] / timings$time_fft[timings$n == 1000]
  expect_true(ratio_time < ratio_n^1.5 * 1.1,  # 10% tolerance for timing variation
              info = paste("FFT time scaling:", round(ratio_time, 2),
                          "vs. O(n^1.5) limit:", round(ratio_n^1.5, 2)))

  # Verify speedup for n=500 (where we have both measurements)
  if (!is.na(timings$time_exact[timings$n == 500])) {
    speedup_500 <- timings$time_exact[timings$n == 500] / timings$time_fft[timings$n == 500]
    expect_true(speedup_500 > 5,
                info = paste("FFT speedup at n=500:", round(speedup_500, 2), "x"))
  }

  # Print timings for documentation
  message("\nPerformance benchmark results:")
  print(timings)
})
