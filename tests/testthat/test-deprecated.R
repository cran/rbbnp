test_that("plot_ft shows deprecation warning", {
  set.seed(123)
  X <- rnorm(50)
  xi_interval <- list(xi_lb = 1, xi_ub = 10)

  # Should warn about deprecation
  expect_warning(
    plot_ft(X, xi_interval = xi_interval),
    "plot_ft.*deprecated"
  )
})


test_that("plot_ft still works after deprecation warning", {
  set.seed(123)
  X <- rnorm(50)
  xi_interval <- list(xi_lb = 1, xi_ub = 10)

  # Should still produce a plot despite warning
  suppressWarnings({
    p <- plot_ft(X, xi_interval = xi_interval)
  })

  expect_s3_class(p, "ggplot")
})
