#' Plot the Fourier Transform
#'
#' Plot the Fourier Transform of the
#'
#' @param X A numerical vector of sample data.
#' @param xi_interval A list containing the lower (`xi_lb`) and upper (`xi_ub`) bounds of the xi interval.
#' @param ft_plot.resol An integer representing the resolution of the plot, specifically the number of points
#'        used to represent the Fourier transform. Defaults to 500.
#'
#' @return A ggplot object representing the plot of the Fourier transform.
#'
#' @details
#' C = 1, the parameter in \eqn{O(1/n^{0.25})}, see more details in Schennach (2020) \doi{10.1093/restud/rdz065}.
#'
#' @export
#'
#' @examples
#' plot_ft(
#'   sample_data$X,
#'   xi_interval = list(xi_lb = 1, xi_ub = 50),
#'   ft_plot.resol = 1000
#' )
plot_ft <- function(X,
                    xi_interval,
                    ft_plot.resol = 500) {

  xi_lb <- xi_interval$xi_lb
  xi_ub <- xi_interval$xi_ub

  xi_ci <- seq(xi_lb, xi_ub, length.out = log(xi_ub / xi_lb) * ft_plot.resol)

  avg_phi <- purrr::map_dbl(xi_ci, function(u) get_avg_phi(xi = u, X = X))

  log_xi = log(xi_ci)
  log_avg_phi = log(avg_phi)
  dt <- data.frame(log_xi, log_avg_phi)
  ggplot(dt, aes(x = log_xi, y = log_avg_phi)) +
    geom_line() +
    geom_vline(xintercept = c(log(xi_lb), log(xi_ub)), color = "grey") +
    annotate("text", x = max(dt$log_xi) * 0.9, y = max(dt$log_avg_phi), label = paste("N =", length(X)), size = 4)
}

#' Bias bound approach for density estimation
#'
#' Estimates the density at a given point or across a range, and provides visualization options for density,
#' bias, and confidence intervals.
#'
#' @param X A numerical vector of sample data.
#' @param x Optional. A scalar or range of points where the density is estimated. If NULL, a range is automatically generated.
#' @param h A scalar bandwidth parameter.
#' @param alpha Confidence level for intervals. Default is 0.05.
#' @param resol Resolution for the estimation range. Default is 100.
#' @param xi_lb Optional. Lower bound for the interval of Fourier Transform frequency xi. Used for determining the range over which A and r is estimated. If NULL, it is automatically determined based on the methods_get_xi.
#' @param xi_ub Optional. Upper bound for the interval of Fourier Transform frequency xi. Similar to xi_lb, it defines the upper range for A and r estimation. If NULL, the upper bound is determined based on the methods_get_xi.
#' @param methods_get_xi A string specifying the method to automatically determine the xi interval if xi_lb and xi_ub are NULL. Options are "Schennach" and "Schennach_loose". If "Schennach" the range is selected based on the Theorem 2 in Schennach2020, if "Schennach_loose", it is defined by the initial interval given in Theorem 2 without selecting the xi_n.
#' @param if_plot_density Logical. If TRUE, plots the density estimation.
#' @param if_plot_ft Logical. If TRUE, plots the Fourier transform.
#' @param ora_Ar Optional list of oracle values for A and r.
#' @param kernel.fun A string specifying the kernel function to be used. Options are "Schennach2004", "sinc", "normal", "epanechnikov".
#' @param if_approx_kernel Logical. If TRUE, uses approximations for the kernel function.
#' @param kernel.resol The resolution for kernel function approximation. See \code{\link{fun_approx}}.
#' @return A list containing various outputs including estimated values, plots, and intervals.
#' @export
#' @examples
#' \donttest{
#' # Example 1: Specifying x for point estimation with manually selected xi range from 1 to 12
#' biasBound_density(
#'   X = sample_data$X,
#'   x = 1,
#'   h = 0.09,
#'   xi_lb = 1,
#'   xi_ub = 12,
#'   if_plot_ft = TRUE,
#'   kernel.fun = "Schennach2004"
#' )
#'
#' # Example 2: Density estimation with manually selected xi range from 1 to 12 xi_lb and xi_ub
#' # biasBound_density(
#' #   X = sample_data$X,
#' #   h = 0.09,
#' #   xi_lb = 1,
#' #   xi_ub = 12,
#' #   if_plot_ft = FALSE,
#' #   kernel.fun = "Schennach2004"
#' # )
#'
#' # Example 3: Density estimation with automatically selected xi range via Theorem 2 in Schennach 2020
#' #  biasBound_density(
#' #   X = sample_data$X,
#' #   h = 0.09,
#' #   methods_get_xi = "Schennach",
#' #   if_plot_ft = TRUE,
#' #   kernel.fun = "Schennach2004"
#' # )
#' }
biasBound_density <- function(X, x = NULL, h = 0.09, alpha = 0.05, resol = 100,
                              xi_lb = NULL, xi_ub = NULL, methods_get_xi = "Schennach",
                              if_plot_density = TRUE, if_plot_ft = FALSE, ora_Ar = NULL,
                              kernel.fun = "Schennach2004", if_approx_kernel = TRUE, kernel.resol = 1000) {
  # determine the kernel and kernel FT function
  if (kernel.fun == "Schennach2004") {
    # if approx the kernel in Schennach2004
    if (if_approx_kernel) {
      # determine the range of the interpolation
      u_lb <- -(max(X) - min(X)) / h * 10
      u_ub <- (max(X) - min(X)) / h * 10

      inf_k <- function(u) fun_approx(u, u_lb = u_lb, u_ub = u_ub, resol = kernel.resol)
    } else {
      inf_k <- W_kernel
    }
    inf_k_ft <- W_kernel_ft
  }

  if (kernel.fun == "sinc") {
    inf_k <- sinc
    inf_k_ft <- sinc_ft
  }

  if (kernel.fun == "normal") {
    inf_k <- normal_kernel
    inf_k_ft <- normal_kernel_ft
  }

  if (kernel.fun == "epanechnikov") {
    inf_k <- epanechnikov_kernel
    inf_k_ft <- epanechnikov_kernel_ft
  }
  #-------------------------------------------

  # get the interval of xi if it is not specified
  if (is.null(xi_lb) | is.null(xi_ub)) {
    xi_interval <- get_xi_interval(X = X, methods = methods_get_xi)
  } else {
    xi_interval <- list(xi_lb = xi_lb, xi_ub = xi_ub)
  }

  # get the estimation of A, r, bias_1x
  est_Ar <- get_est_Ar(X = X, xi_interval = xi_interval)
  b1x <- get_est_b1x(X = X, h = h, est_Ar = est_Ar, inf_k_ft = inf_k_ft)

  return_list <- list(est_Ar = est_Ar, b1x = b1x)

  # plot the Fourier transformation when estimate the A and r
  if (if_plot_ft) {
    p <- plot_ft(X, xi_interval = xi_interval) + geom_abline(intercept = log(est_Ar[1]), slope = -est_Ar[2], color = "red")
    return_list[["ft_plot"]] <- p
  }

  # if not specify the x, then we create the whole interval and plot the estimation
  if (is.null(x)) {
    x_range <- seq(min(X), max(X), length.out = resol)
    f1x <- purrr::map(x_range, get_avg_f1x, X = X, h = h, inf_k = inf_k) %>% unlist()
    sigma <- purrr::map(x_range, get_sigma, X = X, h = h, inf_k = inf_k) %>% unlist()

    lb_bias <- pmax(f1x - b1x, 0)
    ub_bias <- pmax(f1x + b1x, 0)

    lb <- pmax(f1x - sigma * qnorm(1 - alpha / 2) - b1x, 0)
    ub <- pmax(f1x + sigma * qnorm(1 - alpha / 2) + b1x, 0)

    # plot the density estimation with the bands of bias and sd
    if (if_plot_density) {
      data <- data.frame(X = x_range, f1x = f1x, ub_bias = ub_bias, lb_bias = lb_bias, ub = ub, lb = lb)

      gg <- ggplot(data, aes(x = X, y = f1x)) +
        geom_ribbon(aes(ymin = lb, ymax = ub),
          fill = "green", alpha = 0.5
        ) +
        geom_ribbon(aes(ymin = lb_bias, ymax = ub_bias),
          fill = "orange", alpha = 0.5
        ) +
        geom_line(color = "blue") +
        labs(title = "estimated f(x) and CI", x = "X", y = "f(x)")

      # add oracle bias band
      if (!is.null(ora_Ar)) {
        ora_b1x <- get_est_b1x(h = h, est_Ar = ora_Ar)
        lb_ora_bias <- pmax(f1x - ora_b1x, 0)
        ub_ora_bias <- pmax(f1x + ora_b1x, 0)
        gg <- gg + geom_ribbon(aes(ymin = lb_ora_bias, ymax = ub_ora_bias),
          fill = "red", alpha = 0.5
        )
        return_list[["ora_b1x"]] <- ora_b1x
      }

      return_list[["density_plot"]] <- gg
      return_list[["f1x"]] <- f1x
      return_list[["x_range"]] <- x_range
      return_list[["ub"]] <- ub
      return_list[["lb"]] <- lb
      return_list[["ub_bias"]] <- ub_bias
      return_list[["lb_bias"]] <- lb_bias
      return_list[["sigma"]] <- sigma
    }
  } else { # if specific x is provided then give the point estimation
    f1x <- get_avg_f1x(X, x, h, inf_k = inf_k)
    sigma <- get_sigma(X, x, h, inf_k = inf_k)

    lb <- max(c(f1x - sigma * qnorm(1 - alpha / 2) - b1x, 0))
    ub <- max(c(f1x + sigma * qnorm(1 - alpha / 2) + b1x, 0))

    return_list <- c(return_list, list(f1x = f1x, CI = c(lb = lb, ub = ub)))
  }

  return(return_list)
}



#' Bias bound approach for conditional expectation estimation
#'
#' Estimates the density at a given point or across a range, and provides visualization options for density,
#' bias, and confidence intervals.
#'
#' @param Y A numerical vector of sample data.
#' @param X A numerical vector of sample data.
#' @param x Optional. A scalar or range of points where the density is estimated. If NULL, a range is automatically generated.
#' @param h A scalar bandwidth parameter.
#' @param alpha Confidence level for intervals. Default is 0.05.
#' @param est_Ar Optional list of estimates for A and r. If NULL, they are computed using `get_est_Ar()`.
#' @param resol Resolution for the estimation range. Default is 100.
#' @param xi_lb Optional. Lower bound for the interval of Fourier Transform frequency xi. Used for determining the range over which A and r is estimated. If NULL, it is automatically determined based on the methods_get_xi.
#' @param xi_ub Optional. Upper bound for the interval of Fourier Transform frequency xi. Similar to xi_lb, it defines the upper range for A and r estimation. If NULL, the upper bound is determined based on the methods_get_xi.
#' @param methods_get_xi A string specifying the method to automatically determine the xi interval if xi_lb and xi_ub are NULL. Options are "Schennach" and "Schennach_loose". If "Schennach" the range is selected based on the Theorem 2 in Schennach2020, if "Schennach_loose", it is defined by the initial interval given in Theorem 2 without selecting the xi_n.
#' @param if_plot_conditional_mean Logical. If TRUE, plots the conditional mean estimation.
#' @param if_plot_ft Logical. If TRUE, plots the Fourier transform.
#' @param ora_Ar Optional list of oracle values for A and r.
#' @param kernel.fun A string specifying the kernel function to be used. Options are "Schennach2004", "sinc", "normal", "epanechnikov".
#' @param if_approx_kernel Logical. If TRUE, uses approximations for the kernel function.
#' @param kernel.resol The resolution for kernel function approximation. See \code{\link{fun_approx}}.
#' @return A list containing various outputs including estimated values, plots, and intervals.
#' @export
#' @examples
#' \donttest{
#' # Example 1: point estimation of conditional expectation of Y on X
#' biasBound_condExpectation(
#'  Y = sample_data$Y,
#'  X = sample_data$X,
#'  x = 1,
#'  h = 0.09,
#'  kernel.fun = "Schennach2004"
#' )
#'
#' # Example 2: conditional expectation of Y on X with manually selected range of xi
#' # biasBound_condExpectation(
#' # Y = sample_data$Y,
#' #  X = sample_data$X,
#' #  h = 0.09,
#' #  xi_lb = 1,
#' #  xi_ub = 12,
#' #  kernel.fun = "Schennach2004"
#' # )
#' }
biasBound_condExpectation <- function(Y, X, x = NULL, h = 0.09, alpha = 0.05, est_Ar = NULL, resol = 100,
                                      xi_lb = NULL, xi_ub = NULL, methods_get_xi = "Schennach",
                                      if_plot_ft = FALSE, ora_Ar = NULL, if_plot_conditional_mean = TRUE,
                                      kernel.fun = "Schennach2004", if_approx_kernel = TRUE, kernel.resol = 1000) {
  # regularization of Y and X data structure
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length!")
  }

  # determine the kernel and kernel FT function
  if (kernel.fun == "Schennach2004") {
    # if approx the kernel in Schennach2004
    if (if_approx_kernel) {
      # determine the range of the interpolation
      u_lb <- -(max(X) - min(X)) / h * 10
      u_ub <- (max(X) - min(X)) / h * 10

      inf_k <- function(u) fun_approx(u, u_lb = u_lb, u_ub = u_ub, resol = kernel.resol)
    } else {
      inf_k <- W_kernel
    }
    inf_k_ft <- W_kernel_ft
  }

  if (kernel.fun == "sinc") {
    inf_k <- sinc
    inf_k_ft <- sinc_ft
  }

  if (kernel.fun == "normal") {
    inf_k <- normal_kernel
    inf_k_ft <- normal_kernel_ft
  }

  if (kernel.fun == "epanechnikov") {
    inf_k <- epanechnikov_kernel
    inf_k_ft <- epanechnikov_kernel_ft
  }
  #-------------------------------------------

  # get the interval of xi if it is not specified
  if (is.null(xi_lb) | is.null(xi_ub)) {
    xi_interval <- get_xi_interval(Y = Y, X = X, methods = methods_get_xi)
  } else {
    xi_interval <- list(xi_lb = xi_lb, xi_ub = xi_ub)
  }

  # get the estimation of A, r, B, bias_1x, bias_yx
  est_Ar <- get_est_Ar(X = X, xi_interval = xi_interval)
  b1x <- get_est_b1x(X = X, h = h, est_Ar = est_Ar, inf_k_ft = inf_k_ft)
  est_B <- get_est_B(Y = Y)
  byx <- get_est_byx(Y = Y, X = X, h = h, est_Ar = est_Ar, est_B = est_B, inf_k_ft = inf_k_ft)
  return_list <- list(est_Ar = est_Ar, est_B = est_B, b1x = b1x, byx = byx)

  # plot the Fourier transformation when estimate the A and r
  if (if_plot_ft) {
    p <- plot_ft(X) + geom_abline(intercept = log(est_Ar[1]), slope = -est_Ar[2], color = "red")
    return_list[["ft_plot"]] <- p
  }

  # if not specify the x, then wen create the whole interval and plot the estimation
  if (is.null(x)) {
    x_range <- seq(min(X), max(X), length.out = resol)
    f1x <- purrr::map(x_range, get_avg_f1x, X = X, h = h, inf_k = inf_k) %>% unlist()
    fyx <- purrr::map(x_range, get_avg_fyx, Y = Y, X = X, h = h, inf_k = inf_k) %>% unlist()
    sigma <- purrr::map(x_range, get_sigma, X = X, h = h, inf_k = inf_k) %>% unlist()
    sigma_yx <- purrr::map(x_range, get_sigma_yx, Y = Y, X = X, h = h, inf_k = inf_k) %>% unlist()

    lb <- (fyx - byx) / pmax(f1x + sign(fyx - byx) * b1x, 0) - sigma_yx * qnorm(1 - alpha / 2)
    ub <- (fyx + byx) / pmax(f1x - sign(fyx + byx) * b1x, 0) + sigma_yx * qnorm(1 - alpha / 2)

    conditional_mean_yx <- fyx / f1x

    # plot the conditional mean estitmation and its confidence interval
    if (if_plot_conditional_mean) {
      data <- data.frame(X = x_range, conditional_mean_yx = conditional_mean_yx, Y = Y, ub = ub, lb = lb)

      # Set y-axis limits
      ylower <- min(Y, na.rm = TRUE)
      yupper <- max(Y, na.rm = TRUE)
      ylower <- ylower - 0.05 * (yupper - ylower)
      yupper <- yupper + 0.05 * (yupper - ylower)

      # plot the estimation
      gg <- ggplot() +

        # For the ribbon and line
        geom_ribbon(data = data, aes(x = X, ymin = lb, ymax = ub), fill = "grey", alpha = 0.8) +
        geom_line(data = data, aes(x = X, y = conditional_mean_yx), color = "blue") +

        # For the scatter plot
        geom_point(aes(x = X, y = Y), color = "black", alpha = 0.3) +
        labs(title = "estimated E(Y|X = x) and CI", x = "X", y = "E(Y|X = x)") +
        coord_cartesian(ylim = c(ylower, yupper))

      return_list[["density_plot"]] <- gg
      return_list[["conditional_mean_yx"]] <- conditional_mean_yx
      return_list[["ub"]] <- ub
      return_list[["lb"]] <- lb
      return_list[["sigma"]] <- sigma
      return_list[["sigma_yx"]] <- sigma_yx
    }
  } else { # if specific x is provided then give the point estimation
    f1x <- get_avg_f1x(X, x, h, inf_k = inf_k)
    fyx <- get_avg_fyx(Y = Y, X = X, x = x, h = h, inf_k = inf_k)
    sigma <- get_sigma(X, x, h, inf_k = inf_k)
    sigma_yx <- get_sigma_yx(Y = Y, X = X, x = x, h = h, inf_k = inf_k)
    conditional_mean_yx <- fyx / f1x

    lb <- (fyx - byx) / max(c(f1x + sign(fyx - byx) * b1x, 0)) - sigma_yx * qnorm(1 - alpha / 2)
    ub <- (fyx + byx) / max(c(f1x - sign(fyx + byx) * b1x, 0)) + sigma_yx * qnorm(1 - alpha / 2)

    return_list <- c(return_list, list(
      conditional_mean_yx = conditional_mean_yx,
      f1x = f1x, fyx = fyx, sigma = sigma, sigma_yx = sigma_yx,
      CI = c(lb = lb, ub = ub)
    ))
  }

  return(return_list)
}
