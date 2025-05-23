% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_sample.R
\name{gen_sample_data}
\alias{gen_sample_data}
\title{Generate Sample Data}
\usage{
gen_sample_data(size, dgp, seed = NULL)
}
\arguments{
\item{size}{control the sample size.}

\item{dgp}{data generating process, have options "normal", "chisq", "mixed", "poly", "2_fold_uniform".}

\item{seed}{random seed number.}
}
\value{
A numeric vector of length \code{size}. The elements of the vector
are generated according to the specified \code{dgp}:
\describe{
\item{normal}{Normally distributed values with mean 0 and standard deviation 2.}
\item{chisq}{Chi-squared distributed values with df = 10.}
\item{mixed}{Half normally distributed (mean 0, sd = 2) and half chi-squared distributed (df = 10) values.}
\item{poly}{Values from a polynomial cumulative distribution function on \code{[0,1]}.}
\item{2_fold_uniform}{Sum of two uniformly distributed random numbers.}
}
}
\description{
This function used for generate some sample data for experiment
}
