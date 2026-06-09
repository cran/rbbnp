test_that("kernel functions have expected properties", {
  # Test points
  x <- seq(-2, 2, by = 0.5)

  # Test sinc kernel
  sinc_values <- sinc(x)
  expect_equal(sinc(0), 1)  # Special case at x=0 should be 1
  expect_true(all(abs(sinc_values[x != 0]) <= 1))  # Maximum value for x ≠ 0

  # Test normal kernel
  normal_values <- normal_kernel(x)
  expect_equal(normal_kernel(0), 1/sqrt(2*pi))  # Peak at 0
  expect_true(all(normal_values <= 1/sqrt(2*pi)))  # Maximum value

  # Test Epanechnikov kernel
  epanechnikov_values <- epanechnikov_kernel(x)
  expect_equal(epanechnikov_kernel(0), 0.75)  # Peak at 0
  expect_true(all(epanechnikov_values <= 0.75))  # Maximum value
  expect_true(all(epanechnikov_values[x < -1 | x > 1] == 0))  # Outside [-1, 1] should be 0

  # Check that integrals over [-100, 100] are approximately 1
  # Use stats::integrate() (base R) to avoid relying on extra packages.
  expect_equal(stats::integrate(sinc, -100, 100)$value, 1, tolerance = 1e-2)
  expect_equal(stats::integrate(normal_kernel, -100, 100)$value, 1, tolerance = 1e-5)
  expect_equal(stats::integrate(epanechnikov_kernel, -100, 100)$value, 1, tolerance = 1e-5)
  expect_equal(stats::integrate(W_kernel, -100, 100)$value, 1, tolerance = 1e-5)
})
