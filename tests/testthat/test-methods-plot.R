test_that("plot.bbnp_density creates density plot", {
  # Create a simple bbnp_density object for testing
  set.seed(123)
  X <- rnorm(50)
  x_range <- seq(-2, 2, length.out = 20)

  obj <- new_bbnp_density(
    density = abs(rnorm(20, mean = 0.3, sd = 0.05)),
    x = x_range,
    conf_int = list(
      lower = abs(rnorm(20, mean = 0.2, sd = 0.05)),
      upper = abs(rnorm(20, mean = 0.4, sd = 0.05)),
      conf_level = 0.05
    ),
    bias_bound = list(
      b1x = 0.05,
      est_A = 1.2,
      est_r = 0.5,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = X)
  )

  # Test density plot creation
  p <- plot(obj, type = "density")
  expect_s3_class(p, "ggplot")

  # Test customization parameters
  p2 <- plot(obj, type = "density", fill_ci = "blue", alpha_ci = 0.3)
  expect_s3_class(p2, "ggplot")
})


test_that("plot.bbnp_density creates FT plot", {
  # Create a simple bbnp_density object for testing
  set.seed(123)
  X <- rnorm(50)

  obj <- new_bbnp_density(
    density = abs(rnorm(20, mean = 0.3, sd = 0.05)),
    x = seq(-2, 2, length.out = 20),
    bias_bound = list(
      b1x = 0.05,
      est_A = 1.2,
      est_r = 0.5,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = X)
  )

  # Test FT plot creation
  p <- plot(obj, type = "ft", ft_resol = 100)
  expect_s3_class(p, "ggplot")
})


test_that("plot.bbnp_density errors for point estimation", {
  # Create point estimation object
  obj <- new_bbnp_density(
    estimate = 0.3,
    x = 1.0,
    bias_bound = list(
      b1x = 0.05,
      est_A = 1.2,
      est_r = 0.5
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = rnorm(50))
  )

  # Should error for density plot type
  expect_error(
    plot(obj, type = "density"),
    "Cannot create density plot for point estimation"
  )
})


test_that("plot.bbnp_regression creates regression plot", {
  # Create a simple bbnp_regression object for testing
  set.seed(123)
  X <- rnorm(50)
  Y <- X^2 + rnorm(50, sd = 0.5)
  x_range <- seq(-2, 2, length.out = 20)

  obj <- new_bbnp_regression(
    fitted_values = x_range^2,
    x = x_range,
    conf_int = list(
      lower = x_range^2 - 0.5,
      upper = x_range^2 + 0.5,
      conf_level = 0.05
    ),
    bias_bound = list(
      b1x = 0.05,
      byx = 0.08,
      est_A = 1.2,
      est_r = 0.5,
      est_B = 2.0,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = X, Y = Y)
  )

  # Test regression plot creation
  p <- plot(obj, type = "regression")
  expect_s3_class(p, "ggplot")

  # Test customization parameters
  p2 <- plot(obj, type = "regression", fill_ci = "blue", alpha_ci = 0.5)
  expect_s3_class(p2, "ggplot")
})


test_that("plot.bbnp_regression creates FT plot", {
  # Create a simple bbnp_regression object for testing
  set.seed(123)
  X <- rnorm(50)
  Y <- X^2 + rnorm(50, sd = 0.5)

  obj <- new_bbnp_regression(
    fitted_values = rnorm(20),
    x = seq(-2, 2, length.out = 20),
    bias_bound = list(
      b1x = 0.05,
      byx = 0.08,
      est_A = 1.2,
      est_r = 0.5,
      est_B = 2.0,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = X, Y = Y)
  )

  # Test FT plot creation
  p <- plot(obj, type = "ft", ft_resol = 100)
  expect_s3_class(p, "ggplot")
})


test_that("plot.bbnp_regression errors for point estimation", {
  # Create point estimation object
  obj <- new_bbnp_regression(
    estimate = 1.5,
    x = 1.0,
    bias_bound = list(
      b1x = 0.05,
      byx = 0.08,
      est_A = 1.2,
      est_r = 0.5,
      est_B = 2.0
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = rnorm(50), Y = rnorm(50))
  )

  # Should error for regression plot type
  expect_error(
    plot(obj, type = "regression"),
    "Cannot create regression plot for point estimation"
  )
})


test_that("plot methods handle invalid type argument", {
  set.seed(123)
  X <- rnorm(50)

  obj <- new_bbnp_density(
    density = abs(rnorm(20, mean = 0.3, sd = 0.05)),
    x = seq(-2, 2, length.out = 20),
    bias_bound = list(
      b1x = 0.05,
      est_A = 1.2,
      est_r = 0.5,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = X)
  )

  # Should error for invalid type
  expect_error(
    plot(obj, type = "invalid"),
    "'arg' should be one of"
  )
})


test_that("FT plot errors without xi_interval", {
  set.seed(123)
  X <- rnorm(50)

  obj <- new_bbnp_density(
    density = abs(rnorm(20, mean = 0.3, sd = 0.05)),
    x = seq(-2, 2, length.out = 20),
    bias_bound = list(
      b1x = 0.05,
      est_A = 1.2,
      est_r = 0.5,
      xi_interval = NULL  # No xi_interval
    ),
    bandwidth = 0.1,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = X)
  )

  # Should error without xi_interval
  expect_error(
    plot(obj, type = "ft"),
    "No xi_interval found"
  )
})
