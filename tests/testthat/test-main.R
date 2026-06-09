library(testthat)

test_that("biasBound_density works with basic input", {
  # Generate sample data
  set.seed(123)
  X <- gen_sample_data(100, "normal")

  # Test basic functionality
  result <- biasBound_density(X = X, h = 0.3)

  # Check S3 class
  expect_s3_class(result, "bbnp_density")
  expect_s3_class(result, "bbnp")

  # Check structure - should have list structure with S3 class
  expect_type(result, "list")
  expect_true(all(c("bias_bound", "bandwidth", "n", "kernel", "data") %in% names(result)))

  # Check values are reasonable
  expect_true(result$bias_bound$est_A > 0)  # A should be positive
  expect_true(is.numeric(result$bias_bound$b1x))
  expect_true(result$bandwidth > 0)
})

test_that("biasBound_condExpectation works with basic input", {
  # Generate sample data
  set.seed(123)
  X <- gen_sample_data(100, "normal")
  Y <- -X^2 + 3*X + 10 + rnorm(100)

  # Test basic functionality
  result <- biasBound_condExpectation(Y = Y, X = X, h = 0.2)

  # Check S3 class
  expect_s3_class(result, "bbnp_regression")
  expect_s3_class(result, "bbnp")

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("bias_bound", "bandwidth", "n", "kernel", "data") %in% names(result)))

  # Check values are reasonable
  expect_true(result$bias_bound$est_A > 0)
  expect_true(result$bias_bound$est_B > 0)
  expect_true(is.numeric(result$bias_bound$b1x))
  expect_true(is.numeric(result$bias_bound$byx))
})

test_that("biasBound_density works with automatic bandwidth selection", {
  # Generate sample data
  set.seed(123)
  X <- gen_sample_data(100, "normal")

  # Test with Silverman's rule
  result_silver <- biasBound_density(
    X = X,
    h = NULL,
    h_method = "silverman"
  )

  # Check that bandwidth was selected
  expect_true(is.numeric(result_silver$bandwidth))
  expect_true(result_silver$bandwidth > 0)

  # Test with cross-validation
  # Using a small sample to keep tests fast
  X_small <- X[1:30]
  result_cv <- biasBound_density(
    X = X_small,
    h = NULL,
    h_method = "cv"
  )

  # Check that bandwidth was selected
  expect_true(is.numeric(result_cv$bandwidth))
  expect_true(result_cv$bandwidth > 0)
})

test_that("biasBound_condExpectation works with automatic bandwidth selection", {
  # Generate sample data
  set.seed(123)
  X <- gen_sample_data(50, "normal")
  Y <- -X^2 + 3*X + 10 + rnorm(50)

  # Test with Silverman's rule
  result_silver <- biasBound_condExpectation(
    Y = Y,
    X = X,
    h = NULL,
    h_method = "silverman"
  )

  # Check that bandwidth was selected
  expect_true(is.numeric(result_silver$bandwidth))
  expect_true(result_silver$bandwidth > 0)

  # Test with cross-validation (using smaller sample to keep tests fast)
  X_small <- X[1:20]
  Y_small <- Y[1:20]
  result_cv <- biasBound_condExpectation(
    Y = Y_small,
    X = X_small,
    h = NULL,
    h_method = "cv"
  )

  # Check that bandwidth was selected
  expect_true(is.numeric(result_cv$bandwidth))
  expect_true(result_cv$bandwidth > 0)
})
