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
