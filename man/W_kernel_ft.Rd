% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{W_kernel_ft}
\alias{W_kernel_ft}
\title{Define the Fourier transform of a infinite kernel proposed in Schennach 2004}
\usage{
W_kernel_ft(xi, xi_lb = 0.5, xi_ub = 1.5)
}
\arguments{
\item{xi}{A numerical value or vector representing the frequency domain.}

\item{xi_lb}{The lower bound for the frequency domain. Defaults to 0.5.}

\item{xi_ub}{The upper bound for the frequency domain. Defaults to 1.5.}
}
\value{
A numerical value or vector representing the Fourier transform of the infinite
order kernel at the given frequency/frequencies.
}
\description{
Define the Fourier transform of a infinite kernel proposed in Schennach 2004
}
