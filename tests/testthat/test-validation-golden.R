library(testthat)

# Golden tests: lock down a few deterministic outputs so we notice accidental math changes.
# These should be stable across platforms within reasonable tolerances.

# NOTE: these golden values pin the LEGACY window/envelope math (methods_get_xi="Schennach",
# noise_floor="compact", envelope_use_Y=FALSE, integer_r=FALSE). The package default is now
# "snr"/"auto"/TRUE/integer_r=TRUE (rbbnp 1.1.0); validated in test-window-snr.R + the research study.
# We silence warnings so tests focus on numeric stability.

test_that("golden: density fit returns stable coef + first CI rows", {
  set.seed(1)
  X <- gen_sample_data(80, "normal")

  fit <- suppressWarnings(biasBound_density(X, h = 0.25, kernel.fun = "Schennach2004",
                                            methods_get_xi = "Schennach", noise_floor = "compact",
                                            integer_r = FALSE))

  cf <- coef(fit)
  expect_equal(unname(cf["A"]), 0.2658135, tolerance = 1e-6)
  expect_equal(unname(cf["r"]), 1.8261304, tolerance = 1e-6)
  expect_equal(unname(cf["h"]), 0.25, tolerance = 1e-12)

  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(ci[1:5, "lower"], rep(0, 5), tolerance = 1e-12)
  expect_equal(ci[1:5, "upper"], c(0.03327511, 0.03327511, 0.03327511, 0.04336008, 0.05481402), tolerance = 1e-6)
})

test_that("golden: regression fit returns stable coef; CI computed (may be infinite)", {
  set.seed(2)
  X <- gen_sample_data(60, "normal")
  Y <- -X^2 + 3 * X + rnorm(60)

  fit <- suppressWarnings(biasBound_condExpectation(Y, X, h = 0.3, kernel.fun = "Schennach2004",
                                                    methods_get_xi = "Schennach", noise_floor = "compact",
                                                    envelope_use_Y = FALSE, integer_r = FALSE))

  cf <- coef(fit)
  expect_equal(unname(cf["A"]), 0.1928259, tolerance = 1e-6)
  expect_equal(unname(cf["r"]), 1.3564736, tolerance = 1e-6)
  expect_equal(unname(cf["B"]), 6.0108077, tolerance = 1e-6)
  expect_equal(unname(cf["h"]), 0.3, tolerance = 1e-12)

  ci <- confint(fit)
  expect_true(is.matrix(ci))
  # We allow -Inf/+Inf at this stage (denominator can hit 0 in tails); but no NaNs.
  expect_false(any(is.nan(ci)))
})
