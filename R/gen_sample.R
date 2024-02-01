#' Generate Sample Data
#'
#' This function used for generate some sample data for experiment
#'
#' @param size control the sample size.
#' @param dgp data generating process, have options "normal", "chisq", "mixed", "poly", "2_fold_uniform".
#' @param seed random seed number.
#' @return A numeric vector of length \code{size}. The elements of the vector
#' are generated according to the specified \code{dgp}:
#' \describe{
#'   \item{normal}{Normally distributed values with mean 0 and standard deviation 2.}
#'   \item{chisq}{Chi-squared distributed values with df = 10.}
#'   \item{mixed}{Half normally distributed (mean 0, sd = 2) and half chi-squared distributed (df = 10) values.}
#'   \item{poly}{Values from a polynomial cumulative distribution function on \code{[0,1]}.}
#'   \item{2_fold_uniform}{Sum of two uniformly distributed random numbers.}
#' }
gen_sample_data <- function(size, dgp, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (dgp == "normal") {
    x <- rnorm(size, sd = 2)
  }

  if (dgp == "chisq") {
    x <- rchisq(size, df = 10)
  }

  # t distribution + chisq distribution
  if (dgp == "mixed") {
    x <- c(rnorm(size / 2, sd = 2), rchisq(size / 2, df = 10))
  }

  # poly cdf distribution on [0,1]
  if (dgp == "poly") {
    x <- rpoly01(n = size, k = 3)
  }

  if (dgp == "2_fold_uniform") {
    x <- runif(n = size) + runif(n = size)
  }

  return(x)
}

#' Generate n samples from the distribution
#' @param n The number of samples to generate.
#' @param k The exponent in the distribution function, defaults to 5.
#'
#' @return A vector of `n` samples from the specified polynomial distribution.
#'
#' CDF: f(x) = (x-1)^k + 1
rpoly01 <- function(n, k = 5) {
  # Inverse CDF function
  inv_poly_cdf <- function(p) {
    sign(p - 1) * abs(p - 1)^(1 / k) + 1
  }

  samples <- sapply(runif(n), inv_poly_cdf)

  return(samples)
}

#' True density of 2-fold uniform distribution
#' @param x A numerical value or vector where the true density function is evaluated.
#'
#' @return The value of the true density of the 2-fold uniform distribution at each point in `x`.
#'
true_density_2fold <- function(x) {
  convolution_density <- -abs(x - 1) + 1
  return(convolution_density)
}

#' Sample Data
#'
#'
#' @format A data frame with 1000 rows and 2 variables:
#' \describe{
#'   \item{X}{Numeric vector, generated from 2 fold uniform distribution.}
#'   \item{Y}{Numeric vector, `Y = -X^2 + 3*X + rnorm(1000)*X`.}
#' }
#' @docType data
"sample_data"
