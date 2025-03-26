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

test_that("get_sigma_yx handles negative fyx values", {
  # Setup test data that will produce negative fyx
  X <- c(1, 2, 3)
  Y <- c(-10, -20, -30)  # Large negative Y values
  x <- 2
  h <- 0.1
  
  # Test sigma_yx with data likely to produce negative fyx
  result <- get_sigma_yx(Y, X, x, h, inf_k = normal_kernel)
  expect_equal(result, 0)
  
  # Test with data likely to produce positive fyx
  Y2 <- c(1, 4, 9)  # Positive Y values
  result2 <- get_sigma_yx(Y2, X, x, h, inf_k = normal_kernel)
  expect_true(result2 >= 0)
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

test_that("get_sigma_yx handles negative variance and fyx", {
  # Setup test data
  X <- c(1, 2, 3)
  Y <- c(1, 1, 1)  # Same Y values should give zero variance
  x <- 2
  h <- 0.1
  
  # Test with zero variance
  result <- get_sigma_yx(Y, X, x, h, inf_k = normal_kernel)
  expect_equal(result, 0)
  
  # Test with negative fyx
  Y2 <- c(-10, -20, -30)
  result2 <- get_sigma_yx(Y2, X, x, h, inf_k = normal_kernel)
  expect_equal(result2, 0)
  
  # Test with normal data
  Y3 <- c(1, 4, 9)
  result3 <- get_sigma_yx(Y3, X, x, h, inf_k = normal_kernel)
  expect_true(result3 >= 0)
})