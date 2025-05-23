% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bandwidths.R
\name{cv_bandwidth}
\alias{cv_bandwidth}
\title{Cross-Validation for Bandwidth Selection}
\usage{
cv_bandwidth(
  X,
  h_grid = NULL,
  kernel_func,
  kernel_type = "normal",
  grid_size = 512
)
}
\arguments{
\item{X}{A numerical vector of sample data.}

\item{h_grid}{A numerical vector of bandwidth values to evaluate. If NULL (default), a grid is
automatically generated based on the range and distribution of the data.}

\item{kernel_func}{The kernel function to use for cross-validation.}

\item{kernel_type}{A string identifying the kernel type, used only for reference bandwidth.}

\item{grid_size}{Number of grid points for evaluation. Default is 512.}
}
\value{
A scalar representing the optimal bandwidth that minimizes the cross-validation score.
}
\description{
Implements least-squares cross-validation for bandwidth selection with any kernel function.
Uses the self-convolution approach for accurate estimation of the integral term.
}
\examples{
# Generate sample data
X <- rnorm(100)
# Get optimal bandwidth using cross-validation with a normal kernel
kernel_functions <- create_kernel_functions("normal")
h_opt <- cv_bandwidth(X, kernel_func = kernel_functions$kernel,
                     kernel_type = kernel_functions$kernel_type)
}
