# Tests for the SNR window rule (methods_get_xi = "snr") and the two correctness
# parameters added alongside it: noise_floor and envelope_use_Y.

test_that("get_xi_interval('snr') returns a valid, non-collapsed window", {
  set.seed(11); X <- runif(400) + runif(400)        # 2-fold uniform
  w <- get_xi_interval(X = X, methods = "snr")
  expect_type(w, "list")
  expect_named(w, c("xi_lb", "xi_ub"))
  expect_true(w$xi_lb > 0)
  expect_true(w$xi_ub > w$xi_lb)
  # minimum-window guard: never collapses to (near) zero length
  expect_gte(log(w$xi_ub / w$xi_lb), log(1.3) - 1e-8)
})

test_that("snr is non-vacuous (no fallback warning) and narrower than the loose window", {
  set.seed(12); X <- runif(500) + runif(500)
  expect_no_warning(get_xi_interval(X = X, methods = "snr"))
  w_snr   <- get_xi_interval(X = X, methods = "snr")
  w_loose <- get_xi_interval(X = X, methods = "Schennach_loose")
  expect_lt(w_snr$xi_ub, w_loose$xi_ub)          # snr truncates below n^{1/4}/sd
  expect_equal(w_snr$xi_lb, w_loose$xi_lb)        # same lower cutoff
})

test_that("snr never collapses even at small n", {
  for (n in c(60, 120, 300)) {
    set.seed(n); X <- runif(n) + runif(n)
    w <- get_xi_interval(X = X, methods = "snr")
    expect_true(is.finite(w$xi_ub) && w$xi_ub > w$xi_lb)
  }
})

test_that("larger tau gives a narrower snr window", {
  set.seed(13); X <- runif(800) + runif(800)
  expect_gt(.snr_window(X = X, tau = 2)$xi_ub, .snr_window(X = X, tau = 5)$xi_ub)
  # tau may be a function of n
  wf <- .snr_window(X = X, tau = function(n) sqrt(2 * log(n)))
  expect_true(wf$xi_ub > wf$xi_lb)
})

test_that("unknown methods_get_xi errors clearly", {
  set.seed(14); X <- rnorm(100)
  expect_error(get_xi_interval(X = X, methods = "bogus"), "Unknown methods_get_xi")
})

test_that("biasBound_density works with methods_get_xi='snr' (finite CIs)", {
  set.seed(15); X <- runif(400) + runif(400)
  fit <- biasBound_density(X, h = 0.1, methods_get_xi = "snr")
  expect_s3_class(fit, "bbnp_density")
  expect_true(all(is.finite(fit$conf_int$lower)))
  expect_true(all(is.finite(fit$conf_int$upper)))
  expect_true(is.finite(fit$bias_bound$est_r) && fit$bias_bound$est_r > 0)
})

test_that("biasBound_condExpectation works with snr + fixes (auto floor, Y-envelope)", {
  set.seed(16); X <- runif(500) + runif(500); Y <- -X^2 + 3 * X + rnorm(500) * X
  fit <- biasBound_condExpectation(Y, X, h = 0.1, methods_get_xi = "snr",
                                   noise_floor = "auto", envelope_use_Y = TRUE)
  expect_s3_class(fit, "bbnp_regression")
  expect_true(is.finite(fit$bias_bound$est_r))
  expect_true(is.finite(fit$bias_bound$byx))
})

test_that("noise_floor: general > compact for non-constant Y; default is compact", {
  set.seed(17); Y <- rnorm(300, mean = 2); n <- length(Y)
  cmp <- .noise_floor_avg_dphi(Y, n, "compact")
  gen <- .noise_floor_avg_dphi(Y, n, "general")
  expect_gt(gen, cmp)                                   # general form is larger
  expect_equal(.noise_floor_avg_dphi(Y, n, "auto"), gen) # auto -> general for non-constant Y
  expect_equal(.noise_floor_avg_dphi(1, n, "auto"),
               .noise_floor_avg_dphi(1, n, "compact"))  # auto -> compact for constant Y
})

test_that("envelope_use_Y: default ignores Y; TRUE uses phi_{Y;X}", {
  set.seed(18); X <- runif(600) + runif(600); Y <- -X^2 + 3 * X + rnorm(600) * X
  w <- get_xi_interval(X = X, methods = "snr")
  ar_density <- get_est_Ar(Y = 1, X = X, xi_interval = w)
  ar_useF    <- get_est_Ar(Y = Y, X = X, xi_interval = w, use_Y = FALSE)
  ar_useT    <- get_est_Ar(Y = Y, X = X, xi_interval = w, use_Y = TRUE)
  expect_equal(ar_useF, ar_density)                    # default: Y ignored (fit to |phi_X|)
  expect_false(isTRUE(all.equal(ar_useT, ar_useF)))    # use_Y: different (fit to |phi_{Y;X}|)
})

test_that("degenerate (below-floor) envelope warns", {
  set.seed(2); X <- gen_sample_data(60, "normal"); Y <- -X^2 + 3 * X + rnorm(60)
  w <- suppressWarnings(get_xi_interval(Y = Y, X = X, methods = "snr"))
  expect_warning(get_est_Ar(Y = Y, X = X, xi_interval = w, use_Y = TRUE), "smoothness floor")
})

test_that("integer_r=TRUE clamps r>=2 and tames the degenerate bias", {
  set.seed(2); X <- gen_sample_data(60, "normal"); Y <- -X^2 + 3 * X + rnorm(60)
  w <- suppressWarnings(get_xi_interval(Y = Y, X = X, methods = "snr"))
  ar_raw <- suppressWarnings(get_est_Ar(Y, X, w, use_Y = TRUE))
  ar_c1  <- suppressWarnings(get_est_Ar(Y, X, w, use_Y = TRUE, integer_r = TRUE))
  expect_lt(ar_raw[["est_r"]], 1)                            # raw is degenerate
  expect_gte(ar_c1[["est_r"]], 2)                            # clamped to >= 2
  expect_equal(ar_c1[["est_r"]], round(ar_c1[["est_r"]]))    # integer
  f0 <- suppressWarnings(biasBound_condExpectation(Y, X, h = 0.3, integer_r = FALSE))
  f1 <- suppressWarnings(biasBound_condExpectation(Y, X, h = 0.3, integer_r = TRUE))
  expect_gt(f0$bias_bound$byx, 100)                          # legacy (unclamped) blows up
  expect_lt(f1$bias_bound$byx, 10)                           # C1 tames it
  expect_identical(formals(biasBound_density)$integer_r, TRUE)   # default ON (rbbnp 1.1.0)
})

test_that("package defaults are the validated snr / auto / Y-envelope", {
  expect_identical(formals(get_xi_interval)$methods, "snr")
  expect_identical(formals(create_biasBound_config)$noise_floor, "auto")
  expect_identical(formals(create_biasBound_config)$envelope_use_Y, TRUE)
  expect_identical(formals(biasBound_density)$methods_get_xi, "snr")
  expect_identical(formals(biasBound_condExpectation)$envelope_use_Y, TRUE)
  # legacy behaviour is still reachable explicitly
  set.seed(1); X <- runif(300) + runif(300)
  expect_warning(get_xi_interval(X = X, methods = "Schennach"), "No feasible")
})
