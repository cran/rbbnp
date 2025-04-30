#' Infinite Kernel Function
#' @param u A numerical value or vector where the sinc function is evaluated.
#'
#' @return The value of the sinc function at each point in `u`.
sinc <- function(u) {
  ifelse(u == 0,
         1,  # When u = 0, limit of sin(u)/(pi*u) is 1
         sin(pi * u) / (pi * u))
}

#' Define the closed form FT of the infinite order kernel sin(x)/(pi*x)
#'
#' @param x A numerical value or vector where the Fourier Transform is evaluated.
#'
#' @return The value of the Fourier Transform of the sinc function at each point in `x`.
#'
sinc_ft <- function(x) {
  (sign(x + pi) - sign(x - pi)) / 2
}

#' Define the Fourier transform of a infinite kernel proposed in Schennach 2004
#'
#' @param xi A numerical value or vector representing the frequency domain.
#' @param xi_lb The lower bound for the frequency domain. Defaults to 0.5.
#' @param xi_ub The upper bound for the frequency domain. Defaults to 1.5.
#'
#' @return A numerical value or vector representing the Fourier transform of the infinite
#'         order kernel at the given frequency/frequencies.
#'
W_kernel_ft <- function(xi, xi_lb = 0.5, xi_ub = 1.5) {
  condition1 <- abs(xi) <= xi_lb
  condition2 <- abs(xi) <= xi_ub & abs(xi) >= xi_lb

  value_when_condition1 <- 1
  value_when_condition2 <- 1 / (1 + exp((xi_ub - xi_lb) * ((xi_ub - abs(xi))^(-1) - (abs(xi) - xi_lb)^(-1))))
  value_otherwise <- 0

  return(ifelse(condition1, value_when_condition1,
                ifelse(condition2, value_when_condition2,
                       value_otherwise
                )
  ))
}

#' Define the inverse Fourier transform function of W
#' @param u A numerical value or vector representing the time or space domain.
#' @param L The limit for numerical integration, defines the range of integration as \eqn{[-L, L]}.
#'          Defaults to 10.
#'
#' @return A numerical value or vector representing the inverse Fourier transform of the infinite
#'         order kernel at the given time or space point(s).
#'
W_kernel <- function(u, L = 10) {
  inv_ft <- sapply(u, function(t) {
    integrand_func <- function(xi) {
      W_kernel_ft(xi) * exp(1i * xi * t)
    }

    pracma::Real(1 / (2 * pi) * pracma::integral(integrand_func, -L, L))
  })

  return(inv_ft)
}


#' Normal Kernel Function
#'
#' @param u A numerical value or vector representing the input to the kernel function.
#'
#' @return Returns the value of the Normal kernel function at the given input.
#'
normal_kernel <- function(u) {
  dnorm(u, mean = 0, sd = 1)
}

#' Fourier Transform of Normal Kernel
#'
#' @param xi A numerical value or vector representing the frequency domain.
#'
#' @return Returns the value of the Fourier transform of the Normal kernel at the given frequency/frequencies.
#'
normal_kernel_ft <- function(xi) {
  dnorm(xi, mean = 0, sd = 1) * (sqrt(2 * pi))
}

#' Epanechnikov Kernel
#'
#' @param u A numerical value or vector representing the input to the kernel function.
#'
#' @return Returns the value of the Epanechnikov kernel function at the given input.
#'
epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}

#' Fourier Transform Epanechnikov Kernel
#'
#' @param xi A numerical value or vector representing the frequency domain.
#'
#' @return Returns the value of the Fourier transform of the Epanechnikov kernel at the given frequency/frequencies.
#'
epanechnikov_kernel_ft <- function(xi) {
  (3 * sin(xi)) / (2 * xi) - (3 * (2 * xi * cos(xi) + (-2 + xi^2) * sin(xi))) / (2 * xi^3)
}

#' Create kernel functions based on configuration
#'
#' @param kernel.fun A string specifying the kernel function to be used.
#' @param if_approx_kernel Logical. If TRUE, uses approximations for the kernel function.
#' @param kernel.resol The resolution for kernel function approximation.
#'
#' @return A list containing kernel function, its Fourier transform, and the kernel type
#'
#' @export
create_kernel_functions <- function(kernel.fun = "Schennach2004",
                                   if_approx_kernel = TRUE,
                                   kernel.resol = 1000
                                   ) {

  # Store the kernel type (same as input string or derived from function)
  kernel_type <- kernel.fun

  if (kernel.fun == "Schennach2004") {
    if (if_approx_kernel) {
      inf_k <- function(u) fun_approx(u, u_lb = -100, u_ub = 100, resol = kernel.resol)
    } else {
      inf_k <- W_kernel
    }
    inf_k_ft <- W_kernel_ft
  } else if (kernel.fun == "sinc") {
    inf_k <- sinc
    inf_k_ft <- sinc_ft
  } else if (kernel.fun == "normal") {
    inf_k <- normal_kernel
    inf_k_ft <- normal_kernel_ft
  } else if (kernel.fun == "epanechnikov") {
    inf_k <- epanechnikov_kernel
    inf_k_ft <- epanechnikov_kernel_ft
  } else {
    stop("Unsupported kernel function: ", kernel.fun)
  }

  return(list(kernel = inf_k, kernel_ft = inf_k_ft, kernel_type = kernel_type))
}
