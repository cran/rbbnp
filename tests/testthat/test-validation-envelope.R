library(testthat)

# Invariant test: when LP method is used, the estimated (A, r) should satisfy the envelope constraints
# on the grid used internally (within small numerical tolerance).

test_that("LP A/r estimate satisfies envelope constraints on the internal grid", {
  skip_if_not_installed("Rglpk")

  set.seed(123)
  X <- gen_sample_data(120, "normal")

  # Build config so we have xi_interval and the same method plumbing
  cfg <- create_biasBound_config(
    X = X,
    h = 0.2,
    h_method = "cv",
    alpha = 0.05,
    resol = 50,
    xi_lb = NULL,
    xi_ub = NULL,
    methods_get_xi = "Schennach",
    kernel.fun = "Schennach2004",
    if_approx_kernel = TRUE,
    kernel.resol = 300
  )

  xi_interval <- cfg$xi_interval
  est <- suppressWarnings(get_est_Ar(method = "lp", X = X, xi_interval = xi_interval))
  A <- unname(est["est_A"])
  r <- unname(est["est_r"])

  ln_xi_lb <- log(xi_interval$xi_lb)
  ln_xi_ub <- log(xi_interval$xi_ub)
  ln_xi_range <- seq(ln_xi_lb, ln_xi_ub, length.out = (ln_xi_ub - ln_xi_lb) * 200)
  avg_phi_log_values <- vapply(ln_xi_range, function(x) get_avg_phi_log(X = X, ln_xi = x), numeric(1))

  # constraint: log(A) - r * log(xi) >= log(|phi(xi)|)
  lhs <- log(A) - r * ln_xi_range
  expect_true(all(lhs >= avg_phi_log_values - 1e-6))
})
