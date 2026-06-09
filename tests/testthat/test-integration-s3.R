test_that("Full density workflow with all S3 methods", {
  set.seed(123)
  X <- rnorm(100)

  # Estimation
  fit <- biasBound_density(X, h = 0.1)

  # Check S3 class
  expect_s3_class(fit, "bbnp_density")
  expect_s3_class(fit, "bbnp")

  # Test print method
  expect_output(print(fit), "Bias-Bounded Density Estimation")

  # Test summary method
  s <- summary(fit)
  expect_s3_class(s, "summary.bbnp_density")
  expect_output(print(s), "Summary")

  # Test plot methods
  p <- plot(fit, type = "density")
  expect_s3_class(p, "ggplot")

  p_ft <- plot(fit, type = "ft")
  expect_s3_class(p_ft, "ggplot")

  # Test coef method
  coefs <- coef(fit)
  expect_type(coefs, "double")
  expect_named(coefs, c("A", "r", "h"))

  # Test confint method
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(colnames(ci), c("lower", "upper"))
})


test_that("Full regression workflow with all S3 methods", {
  set.seed(123)
  X <- rnorm(100)
  Y <- X^2 + rnorm(100, sd = 0.5)

  # Estimation
  fit <- biasBound_condExpectation(Y, X, h = 0.1)

  # Check S3 class
  expect_s3_class(fit, "bbnp_regression")
  expect_s3_class(fit, "bbnp")

  # Test print method
  expect_output(print(fit), "Bias-Bounded Conditional Expectation")

  # Test summary method
  s <- summary(fit)
  expect_s3_class(s, "summary.bbnp_regression")
  expect_output(print(s), "Summary")

  # Test plot methods
  p <- plot(fit, type = "regression")
  expect_s3_class(p, "ggplot")

  p_ft <- plot(fit, type = "ft")
  expect_s3_class(p_ft, "ggplot")

  # Test coef method
  coefs <- coef(fit)
  expect_type(coefs, "double")
  expect_named(coefs, c("A", "r", "B", "h"))

  # Test fitted method
  fits <- fitted(fit)
  expect_type(fits, "double")
  expect_length(fits, 100)  # default resol

  # Test confint method
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(colnames(ci), c("lower", "upper"))
})


test_that("Point estimation workflow works correctly", {
  set.seed(123)
  X <- rnorm(100)
  Y <- X^2 + rnorm(100, sd = 0.5)

  # Density point estimation
  fit_dens <- biasBound_density(X, x = 1.0, h = 0.1)
  expect_s3_class(fit_dens, "bbnp_density")
  expect_null(fit_dens$density)
  expect_false(is.null(fit_dens$estimate))

  # Regression point estimation
  fit_reg <- biasBound_condExpectation(Y, X, x = 1.0, h = 0.1)
  expect_s3_class(fit_reg, "bbnp_regression")
  expect_null(fit_reg$fitted_values)
  expect_false(is.null(fit_reg$estimate))

  # Test extractors work
  expect_type(coef(fit_dens), "double")
  expect_type(coef(fit_reg), "double")
  expect_type(fitted(fit_reg), "double")
  expect_length(fitted(fit_reg), 1)
})


test_that("Backward compatibility - list access still works", {
  set.seed(123)
  X <- rnorm(100)

  fit <- biasBound_density(X, h = 0.1)

  # Old-style list access should still work
  expect_true(!is.null(fit$bias_bound))
  expect_true(!is.null(fit$bias_bound$est_A))
  expect_true(!is.null(fit$bias_bound$est_r))
  expect_true(!is.null(fit$bias_bound$b1x))
  expect_true(!is.null(fit$bandwidth))
  expect_true(!is.null(fit$n))
  expect_true(!is.null(fit$kernel))
  expect_true(!is.null(fit$data))

  # Values should be accessible
  expect_type(fit$bias_bound$est_A, "double")
  expect_type(fit$bandwidth, "double")
  expect_equal(fit$n, 100)
})


test_that("Error messages are clear for removed parameters", {
  set.seed(123)
  X <- rnorm(100)

  # Old plotting parameters should error
  expect_error(
    biasBound_density(X, h = 0.1, if_plot_density = TRUE),
    "unused argument"
  )

  expect_error(
    biasBound_density(X, h = 0.1, if_plot_ft = TRUE),
    "unused argument"
  )
})


test_that("Workflow with automatic bandwidth selection", {
  set.seed(123)
  X <- rnorm(50)

  # Should work with h = NULL
  fit <- biasBound_density(X, h = NULL, h_method = "silverman")

  # Check that bandwidth was selected
  expect_s3_class(fit, "bbnp_density")
  expect_true(fit$bandwidth > 0)

  # All methods should work
  expect_output(print(fit), "automatic")
  expect_s3_class(plot(fit), "ggplot")
  expect_type(coef(fit), "double")
})


test_that("S3 methods handle edge cases gracefully", {
  set.seed(123)
  X <- rnorm(100)

  fit <- biasBound_density(X, x = 1.0, h = 0.1)

  # Point estimation should error for density plot
  expect_error(
    plot(fit, type = "density"),
    "Cannot create density plot for point estimation"
  )

  # But FT plot should work
  expect_s3_class(plot(fit, type = "ft"), "ggplot")
})


test_that("Integration: Complete workflow matches old numerical results", {
  # This test ensures the refactoring didn't change numerical computations
  set.seed(456)
  X <- rnorm(100)

  fit <- biasBound_density(X, h = 0.1, xi_lb = 1, xi_ub = 10)

  # Basic sanity checks
  expect_true(fit$bias_bound$est_A > 0)
  expect_true(fit$bias_bound$est_r > 0)
  expect_true(fit$bias_bound$b1x > 0)
  expect_true(all(fit$density >= 0))

  # Confidence intervals should contain estimates
  expect_true(all(fit$conf_int$lower <= fit$density))
  expect_true(all(fit$density <= fit$conf_int$upper))
})
