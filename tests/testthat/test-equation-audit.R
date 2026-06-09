library(testthat)

# Equation-level audits: recompute quantities from stored components and
# verify they match what's returned. These tests are designed to catch
# silent discrepancies between code and the mathematical formulas in the
# RJournal draft / theory vignette.


test_that("density CI equals equation: fhat ± b ± z*sigma (then truncated at 0)", {
  set.seed(101)
  X <- gen_sample_data(150, "2_fold_uniform")

  fit <- suppressWarnings(biasBound_density(
    X,
    h = 0.10,
    alpha = 0.05,
    kernel.fun = "Schennach2004",
    resol = 80
  ))

  stopifnot(!is.null(fit$density), !is.null(fit$std_error))

  alpha <- fit$conf_int$conf_level
  z <- qnorm(1 - alpha / 2)
  b <- fit$bias_bound$b1x

  fhat <- fit$density
  sigma <- fit$std_error

  lb_calc <- pmax(fhat - z * sigma - b, 0)
  ub_calc <- pmax(fhat + z * sigma + b, 0)

  expect_equal(fit$conf_int$lower, lb_calc, tolerance = 1e-12)
  expect_equal(fit$conf_int$upper, ub_calc, tolerance = 1e-12)
})


test_that("regression CI equals equation in code (ratio + bias bounds + z*sigma)", {
  set.seed(202)
  n <- 220
  X <- gen_sample_data(n, "2_fold_uniform")
  true_m <- function(x) sin(2 * pi * x)
  Y <- true_m(X) + rnorm(n, sd = 0.3)

  fit <- suppressWarnings(biasBound_condExpectation(
    Y, X,
    h = 0.10,
    alpha = 0.05,
    kernel.fun = "Schennach2004",
    resol = 80
  ))

  alpha <- fit$conf_int$conf_level
  z <- qnorm(1 - alpha / 2)

  fyx <- fit$joint_density
  f1x <- fit$marginal_density
  sigma_yx <- fit$std_error

  b1x <- fit$bias_bound$b1x
  byx <- fit$bias_bound$byx

  lb_calc <- (fyx - byx) / pmax(f1x + sign(fyx - byx) * b1x, 0) - sigma_yx * z
  ub_calc <- (fyx + byx) / pmax(f1x - sign(fyx + byx) * b1x, 0) + sigma_yx * z

  expect_equal(fit$conf_int$lower, lb_calc, tolerance = 1e-12)
  expect_equal(fit$conf_int$upper, ub_calc, tolerance = 1e-12)
})


test_that("bias bound integral has correct order of magnitude (trapz approximation)", {
  skip_on_cran()

  set.seed(303)
  X <- gen_sample_data(200, "2_fold_uniform")
  fit <- suppressWarnings(biasBound_density(
    X,
    h = 0.10,
    alpha = 0.05,
    kernel.fun = "Schennach2004",
    resol = 60
  ))

  # Use estimated envelope parameters
  A <- unname(fit$bias_bound$est_A)
  r <- unname(fit$bias_bound$est_r)
  h <- unname(fit$bandwidth)

  # Access kernel Fourier transform used internally
  inf_k_ft <- fit$internals$kernel_functions$kernel_ft

  # Approximate integral on a finite symmetric grid
  # (This is only a coarse check; it should be close in practice.)
  xi_max <- 200
  m <- 4001
  xi <- seq(-xi_max, xi_max, length.out = m)

  bound_phi <- pmin(1, A * abs(xi)^(-r))
  bound_phi[!is.finite(bound_phi)] <- 1  # handle xi=0

  integrand <- abs(1 - inf_k_ft(xi * h)) * bound_phi

  approx <- pracma::trapz(xi, integrand) / (2 * pi)

  # Stored b1x should be of similar magnitude; allow slack since this is coarse.
  b1x <- fit$bias_bound$b1x
  expect_true(is.finite(approx) && approx > 0)
  expect_true(is.finite(b1x) && b1x > 0)

  # Rough consistency (within an order of magnitude)
  ratio <- approx / b1x
  expect_true(ratio > 0.1 && ratio < 10)
})
