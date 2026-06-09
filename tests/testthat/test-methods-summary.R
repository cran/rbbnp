test_that("summary.bbnp_density works for range estimation", {
  obj <- new_bbnp_density(
    density = c(0.1, 0.2, 0.3),
    x = c(1, 2, 3),
    conf_int = list(conf_level = 0.05),
    bias_bound = list(b1x = 0.01, est_A = 1.5, est_r = 2.0, xi_interval = list(xi_lb = 1, xi_ub = 10)),
    call = quote(biasBound_density(X = rnorm(100))),
    bandwidth = 0.1,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100)),
    internals = list()
  )

  summ <- summary(obj)
  expect_s3_class(summ, "summary.bbnp_density")
  expect_equal(summ$n, 100)
  expect_equal(summ$bandwidth, 0.1)
  expect_false(summ$is_point_est)
  expect_equal(names(summ$density_summary), c("min", "Q1.25%", "median", "mean", "Q3.75%", "max"))
})

test_that("summary.bbnp_density works for point estimation", {
  obj <- new_bbnp_density(
    estimate = 0.25,
    x = 1.5,
    conf_int = list(lower = 0.20, upper = 0.30, conf_level = 0.05),
    bias_bound = list(b1x = 0.01, est_A = 1.5, est_r = 2.0),
    call = quote(biasBound_density(X = rnorm(100), x = 1.5)),
    bandwidth = 0.1,
    n = 100,
    kernel = "normal",
    data = list(X = rnorm(100)),
    internals = list()
  )

  summ <- summary(obj)
  expect_s3_class(summ, "summary.bbnp_density")
  expect_true(summ$is_point_est)
  expect_equal(summ$estimate, 0.25)
})

test_that("print.summary.bbnp_density produces output", {
  obj <- new_bbnp_density(
    density = c(0.1, 0.2, 0.3),
    x = c(1, 2, 3),
    conf_int = list(conf_level = 0.05),
    bias_bound = list(b1x = 0.01, est_A = 1.5, est_r = 2.0),
    call = quote(biasBound_density(X = rnorm(100))),
    bandwidth = 0.1,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100)),
    internals = list()
  )

  output <- capture.output(summary(obj))
  expect_true(any(grepl("Summary: Bias-Bounded Density Estimation", output)))
  expect_true(any(grepl("Sample size", output)))
})

test_that("summary.bbnp_regression works for range estimation", {
  obj <- new_bbnp_regression(
    fitted_values = c(1.1, 1.2, 1.3),
    x = c(1, 2, 3),
    conf_int = list(conf_level = 0.05),
    bias_bound = list(b1x = 0.01, byx = 0.02, est_A = 1.5, est_r = 2.0, est_B = 3.0),
    marginal_density = c(0.3, 0.35, 0.32),
    call = quote(biasBound_condExpectation(Y = rnorm(100), X = rnorm(100))),
    bandwidth = 0.1,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100), Y = rnorm(100)),
    internals = list()
  )

  summ <- summary(obj)
  expect_s3_class(summ, "summary.bbnp_regression")
  expect_equal(summ$n, 100)
  expect_false(summ$is_point_est)
  expect_equal(names(summ$fitted_summary), c("min", "Q1.25%", "median", "mean", "Q3.75%", "max"))
})

test_that("summary.bbnp_regression works for point estimation", {
  obj <- new_bbnp_regression(
    estimate = 1.25,
    x = 1.5,
    conf_int = list(lower = 1.0, upper = 1.5, conf_level = 0.05),
    bias_bound = list(b1x = 0.01, byx = 0.02, est_A = 1.5, est_r = 2.0, est_B = 3.0),
    call = quote(biasBound_condExpectation(Y = rnorm(100), X = rnorm(100), x = 1.5)),
    bandwidth = 0.1,
    n = 100,
    kernel = "normal",
    data = list(X = rnorm(100), Y = rnorm(100)),
    internals = list()
  )

  summ <- summary(obj)
  expect_s3_class(summ, "summary.bbnp_regression")
  expect_true(summ$is_point_est)
  expect_equal(summ$estimate, 1.25)
})
