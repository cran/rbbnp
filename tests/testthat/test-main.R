library(testthat)

test_that("biasBound_density works with basic input", {
  # Generate sample data
  set.seed(123)
  X <- gen_sample_data(100, "normal")

  # Test basic functionality
  result <- biasBound_density(X = X, h = 0.3)

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("est_Ar", "b1x") %in% names(result)))

  # Check values are reasonable
  expect_true(result$est_Ar[1] > 0)  # A should be positive
  expect_true(is.numeric(result$b1x))
})

test_that("biasBound_condExpectation works with basic input", {
  # Generate sample data
  set.seed(123)
  X <- gen_sample_data(100, "normal")
  Y <- -X^2 + 3*X + 10 + rnorm(100)

  # Test basic functionality
  result <- biasBound_condExpectation(Y = Y, X = X, h = 0.2)

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("est_Ar", "est_B", "b1x", "byx") %in% names(result)))
})

test_that("biasBound_density works with automatic bandwidth selection", {
  # Generate sample data
  set.seed(123)
  X <- gen_sample_data(100, "normal")
  
  # Test with Silverman's rule
  result_silver <- biasBound_density(
    X = X, 
    h = NULL, 
    h_method = "silverman",
    if_plot_density = FALSE
  )
  
  # Check that bandwidth was selected
  expect_true(is.numeric(result_silver$h))
  expect_true(result_silver$h > 0)
  
  # Test with cross-validation
  # Using a small sample to keep tests fast
  X_small <- X[1:30]
  result_cv <- biasBound_density(
    X = X_small, 
    h = NULL, 
    h_method = "cv",
    if_plot_density = FALSE
  )
  
  # Check that bandwidth was selected
  expect_true(is.numeric(result_cv$h))
  expect_true(result_cv$h > 0)
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
    h_method = "silverman",
    if_plot_conditional_mean = FALSE
  )
  
  # Check that bandwidth was selected
  expect_true(is.numeric(result_silver$h))
  expect_true(result_silver$h > 0)
  
  # Test with cross-validation (using smaller sample to keep tests fast)
  X_small <- X[1:20]
  Y_small <- Y[1:20]
  result_cv <- biasBound_condExpectation(
    Y = Y_small,
    X = X_small, 
    h = NULL, 
    h_method = "cv",
    if_plot_conditional_mean = FALSE
  )
  
  # Check that bandwidth was selected
  expect_true(is.numeric(result_cv$h))
  expect_true(result_cv$h > 0)
})
