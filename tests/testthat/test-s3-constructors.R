test_that("new_bbnp_density creates valid objects", {
  # Valid range estimation object
  obj <- new_bbnp_density(
    density = c(0.1, 0.2, 0.3),
    x = c(1, 2, 3),
    conf_int = list(lower = c(0.05, 0.15, 0.25), upper = c(0.15, 0.25, 0.35), conf_level = 0.05),
    bias_bound = list(b1x = 0.01, est_A = 1.5, est_r = 2.0, xi_interval = list(xi_lb = 1, xi_ub = 10)),
    std_error = c(0.01, 0.02, 0.03),
    call = quote(biasBound_density(X = rnorm(100))),
    bandwidth = 0.1,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100)),
    internals = list(config = list())
  )

  expect_s3_class(obj, "bbnp_density")
  expect_s3_class(obj, "bbnp")
  expect_equal(length(obj$density), 3)
  expect_equal(obj$bandwidth, 0.1)
  expect_equal(obj$n, 100)
})

test_that("new_bbnp_density creates valid point estimation objects", {
  # Valid point estimation object
  obj <- new_bbnp_density(
    estimate = 0.25,
    x = 1.5,
    conf_int = list(lower = 0.20, upper = 0.30, conf_level = 0.05),
    bias_bound = list(b1x = 0.01, est_A = 1.5, est_r = 2.0, xi_interval = list(xi_lb = 1, xi_ub = 10)),
    std_error = 0.02,
    call = quote(biasBound_density(X = rnorm(100), x = 1.5)),
    bandwidth = 0.1,
    n = 100,
    kernel = "normal",
    data = list(X = rnorm(100)),
    internals = list()
  )

  expect_s3_class(obj, "bbnp_density")
  expect_equal(obj$estimate, 0.25)
  expect_equal(obj$x, 1.5)
})

test_that("new_bbnp_density validates inputs", {
  # Non-numeric density
  expect_error(
    new_bbnp_density(density = "not numeric"),
    "density must be numeric"
  )

  # Non-numeric x
  expect_error(
    new_bbnp_density(x = "not numeric"),
    "x must be numeric"
  )

  # Non-scalar estimate
  expect_error(
    new_bbnp_density(estimate = c(0.1, 0.2)),
    "estimate must be a numeric scalar"
  )

  # Non-positive bandwidth
  expect_error(
    new_bbnp_density(bandwidth = -0.1),
    "bandwidth must be a positive numeric scalar"
  )

  # Non-positive n
  expect_error(
    new_bbnp_density(n = 0),
    "n must be a positive integer"
  )
})

test_that("new_bbnp_regression creates valid objects", {
  # Valid range estimation object
  obj <- new_bbnp_regression(
    fitted_values = c(1.1, 1.2, 1.3),
    x = c(1, 2, 3),
    conf_int = list(lower = c(0.9, 1.0, 1.1), upper = c(1.3, 1.4, 1.5), conf_level = 0.05),
    bias_bound = list(
      b1x = 0.01, byx = 0.02,
      est_A = 1.5, est_r = 2.0, est_B = 3.0,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    std_error = c(0.05, 0.06, 0.07),
    marginal_density = c(0.3, 0.35, 0.32),
    joint_density = c(0.33, 0.42, 0.42),
    call = quote(biasBound_condExpectation(Y = rnorm(100), X = rnorm(100))),
    bandwidth = 0.1,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100), Y = rnorm(100)),
    internals = list()
  )

  expect_s3_class(obj, "bbnp_regression")
  expect_s3_class(obj, "bbnp")
  expect_equal(length(obj$fitted_values), 3)
  expect_equal(obj$bandwidth, 0.1)
  expect_equal(obj$n, 100)
})

test_that("new_bbnp_regression creates valid point estimation objects", {
  # Valid point estimation object
  obj <- new_bbnp_regression(
    estimate = 1.25,
    x = 1.5,
    conf_int = list(lower = 1.0, upper = 1.5, conf_level = 0.05),
    bias_bound = list(b1x = 0.01, byx = 0.02, est_A = 1.5, est_r = 2.0, est_B = 3.0),
    std_error = 0.1,
    marginal_density = 0.3,
    joint_density = 0.375,
    call = quote(biasBound_condExpectation(Y = rnorm(100), X = rnorm(100), x = 1.5)),
    bandwidth = 0.1,
    n = 100,
    kernel = "normal",
    data = list(X = rnorm(100), Y = rnorm(100)),
    internals = list()
  )

  expect_s3_class(obj, "bbnp_regression")
  expect_equal(obj$estimate, 1.25)
  expect_equal(obj$x, 1.5)
})

test_that("new_bbnp_regression validates inputs", {
  # Non-numeric fitted_values
  expect_error(
    new_bbnp_regression(fitted_values = "not numeric"),
    "fitted_values must be numeric"
  )

  # Non-numeric marginal_density
  expect_error(
    new_bbnp_regression(marginal_density = "not numeric"),
    "marginal_density must be numeric"
  )

  # Non-numeric joint_density
  expect_error(
    new_bbnp_regression(joint_density = "not numeric"),
    "joint_density must be numeric"
  )

  # Non-scalar estimate
  expect_error(
    new_bbnp_regression(estimate = c(1.1, 1.2)),
    "estimate must be a numeric scalar"
  )

  # Non-positive bandwidth
  expect_error(
    new_bbnp_regression(bandwidth = -0.1),
    "bandwidth must be a positive numeric scalar"
  )

  # Non-positive n
  expect_error(
    new_bbnp_regression(n = -5),
    "n must be a positive integer"
  )
})
