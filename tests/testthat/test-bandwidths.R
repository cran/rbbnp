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
