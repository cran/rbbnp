% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_avg_phi}
\alias{get_avg_phi}
\title{Compute Sample Average of Fourier Transform Magnitude}
\usage{
get_avg_phi(Y = 1, X, xi)
}
\arguments{
\item{Y}{A numerical vector representing the sample data of variable Y.}

\item{X}{A numerical vector representing the sample data of variable X.}

\item{xi}{A single numerical value representing the frequency at which the Fourier transform
is computed.}
}
\value{
Returns the sample estimation of expected Fourier transform at frequency \code{xi}.
}
\description{
Compute Sample Average of Fourier Transform Magnitude
}
