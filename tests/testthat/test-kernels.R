test_that("kernel functions have expected properties", {
  # Test points
  x <- seq(-2, 2, by = 0.5)
  
  # Test sinc kernel
  sinc_values <- sinc(x)
  expect_equal(sinc(0), 1)  # Special case at x=0 should be 1
  expect_true(all(abs(sinc_values[x != 0]) <= 1/pi))  # Maximum value for x ≠ 0
  
  # Test normal kernel
  normal_values <- normal_kernel(x)
  expect_equal(normal_kernel(0), 1)  # Peak at 0
  expect_true(all(normal_values <= 1))  # Maximum value
}) 