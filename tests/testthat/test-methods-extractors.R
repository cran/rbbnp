test_that("coef.bbnp_density extracts correct coefficients", {
  obj <- new_bbnp_density(
    density = rnorm(20),
    x = seq(-2, 2, length.out = 20),
    bias_bound = list(
      b1x = 0.05,
      est_A = 1.234,
      est_r = 0.567
    ),
    bandwidth = 0.123,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = rnorm(50))
  )

  coefs <- coef(obj)

  expect_type(coefs, "double")
  expect_named(coefs, c("A", "r", "h"))
  expect_equal(unname(coefs["A"]), 1.234)
  expect_equal(unname(coefs["r"]), 0.567)
  expect_equal(unname(coefs["h"]), 0.123)
})


test_that("coef.bbnp_regression extracts correct coefficients", {
  obj <- new_bbnp_regression(
    fitted_values = rnorm(20),
    x = seq(-2, 2, length.out = 20),
    bias_bound = list(
      b1x = 0.05,
      byx = 0.08,
      est_A = 1.234,
      est_r = 0.567,
      est_B = 2.345
    ),
    bandwidth = 0.123,
    n = 50,
    kernel = "Schennach2004",
    data = list(X = rnorm(50), Y = rnorm(50))
  )

  coefs <- coef(obj)

  expect_type(coefs, "double")
  expect_named(coefs, c("A", "r", "B", "h"))
  expect_equal(unname(coefs["A"]), 1.234)
  expect_equal(unname(coefs["r"]), 0.567)
  expect_equal(unname(coefs["B"]), 2.345)
  expect_equal(unname(coefs["h"]), 0.123)
})


test_that("fitted.bbnp_regression extracts fitted values for range estimation", {
  fitted_vals <- c(1.1, 1.2, 1.3, 1.4, 1.5)

  obj <- new_bbnp_regression(
    fitted_values = fitted_vals,
    x = seq(-2, 2, length.out = 5),
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

  fit <- fitted(obj)

  expect_type(fit, "double")
  expect_equal(fit, fitted_vals)
  expect_length(fit, 5)
})


test_that("fitted.bbnp_regression extracts estimate for point estimation", {
  obj <- new_bbnp_regression(
    estimate = 1.234,
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

  fit <- fitted(obj)

  expect_type(fit, "double")
  expect_equal(fit, 1.234)
  expect_length(fit, 1)
})


test_that("fitted.bbnp_regression errors when no fitted values exist", {
  obj <- new_bbnp_regression(
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

  expect_error(
    fitted(obj),
    "No fitted values available"
  )
})


test_that("confint.bbnp_density extracts confidence intervals for range estimation", {
  lower_vals <- c(0.1, 0.15, 0.2, 0.25, 0.3)
  upper_vals <- c(0.3, 0.35, 0.4, 0.45, 0.5)

  obj <- new_bbnp_density(
    density = c(0.2, 0.25, 0.3, 0.35, 0.4),
    x = seq(-2, 2, length.out = 5),
    conf_int = list(
      lower = lower_vals,
      upper = upper_vals,
      conf_level = 0.05
    ),
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

  ci <- confint(obj)

  expect_type(ci, "double")
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 5)
  expect_equal(ncol(ci), 2)
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_equal(ci[, "lower"], lower_vals)
  expect_equal(ci[, "upper"], upper_vals)
})


test_that("confint.bbnp_density extracts confidence intervals for point estimation", {
  obj <- new_bbnp_density(
    estimate = 0.3,
    x = 1.0,
    conf_int = list(
      lower = 0.2,
      upper = 0.4,
      conf_level = 0.05
    ),
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

  ci <- confint(obj)

  expect_type(ci, "double")
  expect_false(is.matrix(ci))
  expect_named(ci, c("lower", "upper"))
  expect_equal(unname(ci["lower"]), 0.2)
  expect_equal(unname(ci["upper"]), 0.4)
})


test_that("confint.bbnp_regression extracts confidence intervals for range estimation", {
  lower_vals <- c(1.1, 1.2, 1.3, 1.4, 1.5)
  upper_vals <- c(2.1, 2.2, 2.3, 2.4, 2.5)

  obj <- new_bbnp_regression(
    fitted_values = c(1.6, 1.7, 1.8, 1.9, 2.0),
    x = seq(-2, 2, length.out = 5),
    conf_int = list(
      lower = lower_vals,
      upper = upper_vals,
      conf_level = 0.05
    ),
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

  ci <- confint(obj)

  expect_type(ci, "double")
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 5)
  expect_equal(ncol(ci), 2)
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_equal(ci[, "lower"], lower_vals)
  expect_equal(ci[, "upper"], upper_vals)
})


test_that("confint.bbnp_regression extracts confidence intervals for point estimation", {
  obj <- new_bbnp_regression(
    estimate = 1.5,
    x = 1.0,
    conf_int = list(
      lower = 1.2,
      upper = 1.8,
      conf_level = 0.05
    ),
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

  ci <- confint(obj)

  expect_type(ci, "double")
  expect_false(is.matrix(ci))
  expect_named(ci, c("lower", "upper"))
  expect_equal(unname(ci["lower"]), 1.2)
  expect_equal(unname(ci["upper"]), 1.8)
})


test_that("confint methods error when no CI available", {
  obj_density <- new_bbnp_density(
    density = rnorm(20),
    x = seq(-2, 2, length.out = 20),
    conf_int = list(),
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

  expect_error(
    confint(obj_density),
    "No confidence intervals available"
  )

  obj_regression <- new_bbnp_regression(
    fitted_values = rnorm(20),
    x = seq(-2, 2, length.out = 20),
    conf_int = list(),
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

  expect_error(
    confint(obj_regression),
    "No confidence intervals available"
  )
})


test_that("confint methods warn about ignored level parameter", {
  obj <- new_bbnp_density(
    estimate = 0.3,
    x = 1.0,
    conf_int = list(
      lower = 0.2,
      upper = 0.4,
      conf_level = 0.05
    ),
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

  expect_warning(
    confint(obj, level = 0.90),
    "The 'level' parameter is ignored"
  )
})
