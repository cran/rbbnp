test_that("print.bbnp_density works for range estimation", {
  # Create a sample object for range estimation
  obj <- new_bbnp_density(
    density = c(0.1, 0.2, 0.3),
    x = c(1, 2, 3),
    conf_int = list(
      lower = c(0.05, 0.15, 0.25),
      upper = c(0.15, 0.25, 0.35),
      conf_level = 0.05
    ),
    bias_bound = list(
      b1x = 0.01,
      est_A = 1.5,
      est_r = 2.0,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    std_error = c(0.01, 0.02, 0.03),
    call = quote(biasBound_density(X = rnorm(100), h = 0.1)),
    bandwidth = 0.1,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100)),
    internals = list()
  )

  # Capture output
  output <- capture.output(print(obj))

  # Check key information appears
  expect_true(any(grepl("Bias-Bounded Density Estimation", output)))
  expect_true(any(grepl("Sample size: n = 100", output)))
  expect_true(any(grepl("h = 0.1000", output)))
  expect_true(any(grepl("Kernel:\\s+Schennach2004", output)))
  expect_true(any(grepl("A = 1.5000, r = 2.0000", output)))
  expect_true(any(grepl("b1x = 0.0100", output)))
  expect_true(any(grepl("Evaluation points: 3", output)))
  expect_true(any(grepl("Confidence level: 95%", output)))
})

test_that("print.bbnp_density works for point estimation", {
  # Create a sample object for point estimation
  obj <- new_bbnp_density(
    estimate = 0.25,
    x = 1.5,
    conf_int = list(lower = 0.20, upper = 0.30, conf_level = 0.05),
    bias_bound = list(
      b1x = 0.01,
      est_A = 1.5,
      est_r = 2.0,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    std_error = 0.02,
    call = quote(biasBound_density(X = rnorm(100), x = 1.5, h = 0.1)),
    bandwidth = 0.1,
    n = 100,
    kernel = "normal",
    data = list(X = rnorm(100)),
    internals = list()
  )

  # Capture output
  output <- capture.output(print(obj))

  # Check key information appears
  expect_true(any(grepl("Point estimate at x = 1.5000: f\\(x\\) = 0.2500", output)))
  expect_true(any(grepl("Kernel:\\s+normal", output)))
})

test_that("print.bbnp_density shows automatic bandwidth correctly", {
  # Create object with NULL h in call (automatic selection)
  obj <- new_bbnp_density(
    density = c(0.1, 0.2, 0.3),
    x = c(1, 2, 3),
    conf_int = list(conf_level = 0.05),
    bias_bound = list(b1x = 0.01, est_A = 1.5, est_r = 2.0),
    call = quote(biasBound_density(X = rnorm(100), h = NULL)),
    bandwidth = 0.0987,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100)),
    internals = list()
  )

  output <- capture.output(print(obj))
  expect_true(any(grepl("automatic", output)))
})

test_that("print.bbnp_regression works for range estimation", {
  # Create a sample object for range estimation
  obj <- new_bbnp_regression(
    fitted_values = c(1.1, 1.2, 1.3),
    x = c(1, 2, 3),
    conf_int = list(
      lower = c(0.9, 1.0, 1.1),
      upper = c(1.3, 1.4, 1.5),
      conf_level = 0.05
    ),
    bias_bound = list(
      b1x = 0.01,
      byx = 0.02,
      est_A = 1.5,
      est_r = 2.0,
      est_B = 3.0,
      xi_interval = list(xi_lb = 1, xi_ub = 10)
    ),
    std_error = c(0.05, 0.06, 0.07),
    marginal_density = c(0.3, 0.35, 0.32),
    joint_density = c(0.33, 0.42, 0.42),
    call = quote(biasBound_condExpectation(Y = rnorm(100), X = rnorm(100), h = 0.1)),
    bandwidth = 0.1,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100), Y = rnorm(100)),
    internals = list()
  )

  # Capture output
  output <- capture.output(print(obj))

  # Check key information appears
  expect_true(any(grepl("Bias-Bounded Conditional Expectation Estimation", output)))
  expect_true(any(grepl("Sample size: n = 100", output)))
  expect_true(any(grepl("h = 0.1000", output)))
  expect_true(any(grepl("Kernel:\\s+Schennach2004", output)))
  expect_true(any(grepl("A = 1.5000, r = 2.0000, B = 3.0000", output)))
  expect_true(any(grepl("b1x = 0.0100, byx = 0.0200", output)))
  expect_true(any(grepl("Evaluation points: 3", output)))
  expect_true(any(grepl("Fitted values: E\\[Y\\|X\\]", output)))
  expect_true(any(grepl("Confidence level: 95%", output)))
})

test_that("print.bbnp_regression works for point estimation", {
  # Create a sample object for point estimation
  obj <- new_bbnp_regression(
    estimate = 1.25,
    x = 1.5,
    conf_int = list(lower = 1.0, upper = 1.5, conf_level = 0.05),
    bias_bound = list(
      b1x = 0.01,
      byx = 0.02,
      est_A = 1.5,
      est_r = 2.0,
      est_B = 3.0
    ),
    std_error = 0.1,
    marginal_density = 0.3,
    joint_density = 0.375,
    call = quote(biasBound_condExpectation(Y = rnorm(100), X = rnorm(100), x = 1.5, h = 0.1)),
    bandwidth = 0.1,
    n = 100,
    kernel = "normal",
    data = list(X = rnorm(100), Y = rnorm(100)),
    internals = list()
  )

  # Capture output
  output <- capture.output(print(obj))

  # Check key information appears
  expect_true(any(grepl("Point estimate at x = 1.5000: E\\[Y\\|X=x\\] = 1.2500", output)))
  expect_true(any(grepl("Kernel:\\s+normal", output)))
})

test_that("print.bbnp_regression shows automatic bandwidth correctly", {
  # Create object with NULL h in call (automatic selection)
  obj <- new_bbnp_regression(
    fitted_values = c(1.1, 1.2, 1.3),
    x = c(1, 2, 3),
    conf_int = list(conf_level = 0.05),
    bias_bound = list(b1x = 0.01, byx = 0.02, est_A = 1.5, est_r = 2.0, est_B = 3.0),
    call = quote(biasBound_condExpectation(Y = rnorm(100), X = rnorm(100), h = NULL)),
    bandwidth = 0.0987,
    n = 100,
    kernel = "Schennach2004",
    data = list(X = rnorm(100), Y = rnorm(100)),
    internals = list()
  )

  output <- capture.output(print(obj))
  expect_true(any(grepl("automatic", output)))
})
