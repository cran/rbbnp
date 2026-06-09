# Suppress R CMD check NOTEs for ggplot2 NSE (non-standard evaluation)
# These are column names in data frames used within aes()
utils::globalVariables(c("X", "density", "conditional_mean",
                         "lb", "ub", "lb_bias", "ub_bias",
                         "log_xi", "log_avg_phi",
                         "x", "y", "xend", "yend"))

# Shared palette for the bbnp plot theme (colorblind-friendly).
.bbnp_pal <- c(estimate = "#08306B", ci = "#9ECAE1", bias = "#4292C6",
               point = "grey60", envelope = "#CB181D", window = "#525252")

#' Internal Helper: Minimal bbnp plot theme
#'
#' A clean ggplot2 theme (minimal background, subtle grid, left-aligned bold
#' title) shared by all bbnp plots.
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme object.
#' @keywords internal
#' @noRd
theme_bbnp <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle      = ggplot2::element_text(colour = "grey40", hjust = 0),
      plot.title.position = "plot",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_line(colour = "grey92"),
      axis.title         = ggplot2::element_text(colour = "grey25"),
      axis.text          = ggplot2::element_text(colour = "grey45"),
      legend.position    = "top",
      legend.justification = "left",
      legend.title       = ggplot2::element_blank()
    )
}

#' Plot Method for bbnp_density Objects
#'
#' Creates visualizations of bias-bounded density estimation results
#'
#' @param x An object of class bbnp_density
#' @param type Character string specifying plot type. Options are:
#'   \itemize{
#'     \item \code{"density"} (default): Density estimate with bias bounds and confidence intervals
#'     \item \code{"ft"}: Fourier transform plot with estimated envelope
#'   }
#' @param fill_ci Color for confidence interval ribbon (default: a muted blue).
#' @param fill_bias Color for bias bound ribbon (default: a muted terracotta).
#' @param alpha_ci Transparency for confidence interval ribbon (default: 0.30)
#' @param alpha_bias Transparency for bias bound ribbon (default: 0.45)
#' @param ft_resol Resolution for Fourier transform plot (default: 500)
#' @param xi_range Optional numeric \code{c(lower, upper)} giving the frequency
#'   range to *display* in the \code{"ft"} plot. If \code{NULL} (default) a wide
#'   range around the selected window is shown (see \code{expand}). This controls
#'   only what is drawn; it does not change the fitting window \code{[xi_lb, xi_ub]}.
#' @param expand For the \code{"ft"} plot when \code{xi_range} is \code{NULL}: how
#'   far past the selected window to display, as a multiple (default 2.5).
#' @param ... Additional arguments (unused)
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \donttest{
#' X <- rnorm(100)
#' fit <- biasBound_density(X, h = 0.1)
#' plot(fit)
#' plot(fit, type = "ft")
#' }
plot.bbnp_density <- function(x,
                               type = c("density", "ft"),
                               fill_ci = .bbnp_pal[["ci"]],
                               fill_bias = .bbnp_pal[["bias"]],
                               alpha_ci = 0.55,
                               alpha_bias = 0.55,
                               ft_resol = 500,
                               xi_range = NULL,
                               expand = 2.5,
                               ...) {

  type <- match.arg(type)

  if (type == "density") {
    plot_density_internal(x, fill_ci, fill_bias, alpha_ci, alpha_bias)
  } else if (type == "ft") {
    plot_ft_internal(x, ft_resol, xi_range = xi_range, expand = expand)
  }
}


#' Plot Method for bbnp_regression Objects
#'
#' Creates visualizations of bias-bounded conditional expectation estimation results
#'
#' @param x An object of class bbnp_regression
#' @param type Character string specifying plot type. Options are:
#'   \itemize{
#'     \item \code{"regression"} (default): Conditional expectation with confidence interval
#'     \item \code{"ft"}: Fourier transform plot with estimated envelope
#'   }
#' @param fill_ci Color for confidence interval ribbon (default: a muted blue).
#' @param alpha_ci Transparency for confidence interval ribbon (default: 0.35)
#' @param point_alpha Transparency for data points (default: 0.28)
#' @param point_color Color for data points (default: a soft grey).
#' @param ft_resol Resolution for Fourier transform plot (default: 500)
#' @param xi_range Optional numeric \code{c(lower, upper)} giving the frequency
#'   range to *display* in the \code{"ft"} plot. If \code{NULL} (default) a wide
#'   range around the selected window is shown (see \code{expand}). This controls
#'   only what is drawn; it does not change the fitting window \code{[xi_lb, xi_ub]}.
#' @param expand For the \code{"ft"} plot when \code{xi_range} is \code{NULL}: how
#'   far past the selected window to display, as a multiple (default 2.5).
#' @param ... Additional arguments (unused)
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \donttest{
#' X <- rnorm(100)
#' Y <- X^2 + rnorm(100)
#' fit <- biasBound_condExpectation(Y, X, h = 0.1)
#' plot(fit)
#' plot(fit, type = "ft")
#' }
plot.bbnp_regression <- function(x,
                                  type = c("regression", "ft"),
                                  fill_ci = .bbnp_pal[["ci"]],
                                  alpha_ci = 0.55,
                                  point_alpha = 0.28,
                                  point_color = .bbnp_pal[["point"]],
                                  ft_resol = 500,
                                  xi_range = NULL,
                                  expand = 2.5,
                                  ...) {

  type <- match.arg(type)

  if (type == "regression") {
    plot_regression_internal(x, fill_ci, alpha_ci, point_alpha, point_color)
  } else if (type == "ft") {
    plot_ft_internal(x, ft_resol, xi_range = xi_range, expand = expand)
  }
}


#' Internal Helper: Plot Density Estimation
#'
#' Creates density plot with bias bounds and confidence intervals
#'
#' @param x An object of class bbnp_density
#' @param fill_ci Color for confidence interval ribbon
#' @param fill_bias Color for bias bound ribbon
#' @param alpha_ci Transparency for confidence interval ribbon
#' @param alpha_bias Transparency for bias bound ribbon
#'
#' @return A ggplot2 object
#' @keywords internal
plot_density_internal <- function(x, fill_ci, fill_bias, alpha_ci, alpha_bias) {

  # Check if range estimation exists
  if (is.null(x$density) || is.null(x$x)) {
    stop("Cannot create density plot for point estimation. Use type = 'ft' instead.",
         call. = FALSE)
  }

  # Extract data
  x_vals <- x$x
  density_vals <- x$density

  # Compute bias bounds
  b1x <- x$bias_bound$b1x
  lb_bias <- pmax(density_vals - b1x, 0)
  ub_bias <- pmax(density_vals + b1x, 0)

  # Extract confidence intervals
  lb <- x$conf_int$lower
  ub <- x$conf_int$upper

  # Create data frame for plotting
  data <- data.frame(
    X = x_vals,
    density = density_vals,
    ub_bias = ub_bias,
    lb_bias = lb_bias,
    ub = ub,
    lb = lb
  )

  subtitle <- sprintf("n = %d  |  h = %.3f  |  r = %.0f",
                      x$n, x$bandwidth, x$bias_bound$est_r)

  # Create ggplot
  gg <- ggplot2::ggplot(data, ggplot2::aes(x = X, y = density)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lb, ymax = ub, fill = "95% CI"),
                        alpha = alpha_ci) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lb_bias, ymax = ub_bias, fill = "Bias bound"),
                        alpha = alpha_bias) +
    ggplot2::geom_line(ggplot2::aes(colour = "Estimate"), linewidth = 0.9) +
    ggplot2::scale_fill_manual(values = c("95% CI" = fill_ci, "Bias bound" = fill_bias)) +
    ggplot2::scale_colour_manual(values = c("Estimate" = .bbnp_pal[["estimate"]])) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 2,
                                   override.aes = list(alpha = c(alpha_ci, alpha_bias)))) +
    ggplot2::labs(title = "Bias-bounded density", subtitle = subtitle,
                 x = "x", y = "f(x)") +
    theme_bbnp()

  return(gg)
}


#' Internal Helper: Plot Conditional Expectation
#'
#' Creates regression plot with confidence interval and data points
#'
#' @param x An object of class bbnp_regression
#' @param fill_ci Color for confidence interval ribbon
#' @param alpha_ci Transparency for confidence interval ribbon
#' @param point_alpha Transparency for data points
#' @param point_color Color for data points
#'
#' @return A ggplot2 object
#' @keywords internal
plot_regression_internal <- function(x, fill_ci, alpha_ci, point_alpha, point_color) {

  # Check if range estimation exists
  if (is.null(x$fitted_values) || is.null(x$x)) {
    stop("Cannot create regression plot for point estimation. Use type = 'ft' instead.",
         call. = FALSE)
  }

  # Extract data
  x_vals <- x$x
  fitted_vals <- x$fitted_values

  # Extract confidence intervals
  lb <- x$conf_int$lower
  ub <- x$conf_int$upper

  # Get original data
  X_data <- x$data$X
  Y_data <- x$data$Y

  # Create data frame for plotting
  data <- data.frame(
    X = x_vals,
    conditional_mean = fitted_vals,
    ub = ub,
    lb = lb
  )

  # Set y-axis limits
  ylower <- min(Y_data, na.rm = TRUE)
  yupper <- max(Y_data, na.rm = TRUE)
  ylower <- ylower - 0.05 * (yupper - ylower)
  yupper <- yupper + 0.05 * (yupper - ylower)

  # Where the marginal density is near zero the ratio-based interval can blow up
  # (becoming infinite or far wider than the data). Such a band would simply fill
  # the panel without being informative, so the ribbon is drawn only where the
  # interval is finite and no wider than the visible range.
  span <- yupper - ylower
  ribbon_data <- data[is.finite(data$lb) & is.finite(data$ub) &
                        (data$ub - data$lb) <= span, , drop = FALSE]

  subtitle <- sprintf("n = %d  |  h = %.3f", x$n, x$bandwidth)

  # Create ggplot
  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = X_data, y = Y_data),
                       colour = point_color, alpha = point_alpha,
                       size = 0.9, stroke = 0) +
    ggplot2::geom_ribbon(data = ribbon_data,
                        ggplot2::aes(x = X, ymin = lb, ymax = ub, fill = "95% CI"),
                        alpha = alpha_ci) +
    ggplot2::geom_line(data = data,
                      ggplot2::aes(x = X, y = conditional_mean, colour = "Estimate"),
                      linewidth = 0.9) +
    ggplot2::scale_fill_manual(values = c("95% CI" = fill_ci)) +
    ggplot2::scale_colour_manual(values = c("Estimate" = .bbnp_pal[["estimate"]])) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 2, override.aes = list(alpha = alpha_ci))) +
    ggplot2::labs(title = "Bias-bounded conditional expectation", subtitle = subtitle,
                 x = "x", y = "E[Y | X = x]") +
    ggplot2::coord_cartesian(ylim = c(ylower, yupper)) +
    theme_bbnp()

  return(gg)
}


#' Internal Helper: Plot Fourier Transform
#'
#' Creates a Fourier transform plot over a wide frequency range with the
#' rule-selected fitting window shaded and the estimated envelope overlaid.
#' For regression objects the cross-spectrum `|phi_YX|` is shown (matching the
#' fitted envelope); for density objects the marginal `|phi_X|` is shown.
#'
#' @param x An object of class bbnp_density or bbnp_regression
#' @param ft_resol Resolution for Fourier transform plot
#' @param xi_range Optional display range c(lower, upper); NULL = wide auto range
#' @param expand How far past the window to display when xi_range is NULL
#'
#' @return A ggplot2 object
#' @keywords internal
plot_ft_internal <- function(x, ft_resol, xi_range = NULL, expand = 2.5) {

  # Extract original data
  X <- x$data$X
  n <- length(X)
  sdx <- sqrt(stats::var(X))

  # For regression, use the cross-spectrum |phi_{Y;X}| so the curve matches the
  # fitted envelope; density uses the marginal (Y = 1).
  Y_ft <- if (inherits(x, "bbnp_regression") && !is.null(x$data$Y)) x$data$Y else 1

  # Extract xi interval (the fitting window)
  xi_interval <- x$bias_bound$xi_interval
  if (is.null(xi_interval)) {
    stop("No xi_interval found in object. Cannot create FT plot.",
         call. = FALSE)
  }

  xi_lb <- xi_interval$xi_lb
  xi_ub <- xi_interval$xi_ub

  # Display range: a wide view around the window, or a user override.
  if (is.null(xi_range)) {
    disp_lb <- xi_lb / expand
    disp_ub <- max(n^0.25 / sdx, xi_ub * expand)
  } else {
    disp_lb <- xi_range[1]
    disp_ub <- xi_range[2]
  }

  # Generate xi sequence (log-spaced over the display range)
  xi_ci <- exp(seq(log(disp_lb), log(disp_ub),
                   length.out = max(50, log(disp_ub / disp_lb) * ft_resol)))

  # Compute average phi for each xi
  avg_phi <- purrr::map_dbl(xi_ci, function(u) get_avg_phi(Y = Y_ft, xi = u, X = X))

  # Create data frame
  dt <- data.frame(log_xi = log(xi_ci), log_avg_phi = log(avg_phi))

  # Extract estimated A and r
  est_A <- x$bias_bound$est_A
  est_r <- x$bias_bound$est_r

  subtitle <- sprintf("n = %d  |  window [%.2f, %.2f]  |  A = %.3f, r = %.2f",
                      n, xi_lb, xi_ub, est_A, est_r)

  # Fitted envelope as a segment confined to the selected window
  env_df <- data.frame(
    x = log(xi_lb), xend = log(xi_ub),
    y = log(est_A) - est_r * log(xi_lb),
    yend = log(est_A) - est_r * log(xi_ub)
  )

  # Create ggplot
  gg <- ggplot2::ggplot(dt, ggplot2::aes(x = log_xi, y = log_avg_phi)) +
    ggplot2::annotate("rect", xmin = log(xi_lb), xmax = log(xi_ub),
                     ymin = -Inf, ymax = Inf,
                     fill = .bbnp_pal[["ci"]], alpha = 0.12) +
    ggplot2::geom_vline(xintercept = c(log(xi_lb), log(xi_ub)),
                       colour = .bbnp_pal[["window"]], linetype = 2, linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(colour = "Empirical |phi|"), linewidth = 0.75) +
    ggplot2::geom_segment(data = env_df, inherit.aes = FALSE,
                         ggplot2::aes(x = x, xend = xend, y = y, yend = yend,
                                      colour = "Fitted envelope"), linewidth = 1) +
    ggplot2::annotate("text", x = (log(xi_lb) + log(xi_ub)) / 2, y = Inf,
                     label = "selected window", vjust = 1.5, size = 3.1,
                     colour = .bbnp_pal[["window"]]) +
    ggplot2::scale_colour_manual(values = c("Empirical |phi|" = "grey35",
                                           "Fitted envelope" = .bbnp_pal[["envelope"]])) +
    ggplot2::labs(title = "Fourier transform & fitted envelope", subtitle = subtitle,
                 x = "log |xi|", y = "log |phi(xi)|") +
    theme_bbnp()

  return(gg)
}
