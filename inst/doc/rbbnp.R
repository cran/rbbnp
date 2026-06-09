## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.align = "center",
  out.width = "80%",
  warning = FALSE,
  message = FALSE
)
library(rbbnp)
library(ggplot2)

## ----eval=FALSE---------------------------------------------------------------
# # Install from CRAN
# install.packages("rbbnp")
# 
# # Or install development version from GitHub
# # install.packages("devtools")
# devtools::install_github("xinyu-daidai/rbbnp-dev")

## ----density-quick------------------------------------------------------------
# Generate sample data
X <- gen_sample_data(size = 500, dgp = "2_fold_uniform", seed = 123456)

# Estimate density with bias-aware confidence intervals
fit <- biasBound_density(X, h = 0.1, kernel.fun = "Schennach2004")

# View summary
fit

## ----density-plot-------------------------------------------------------------
# Visualize results
plot(fit)

## ----regression-quick---------------------------------------------------------
# Generate regression data: Y = -X^2 + 3X + noise
Y <- -X^2 + 3*X + rnorm(500) * X

# Estimate conditional expectation
fit_reg <- biasBound_condExpectation(Y, X, h = 0.1, kernel.fun = "Schennach2004")

# Visualize
plot(fit_reg)

## ----methods-demo-------------------------------------------------------------
# Extract parameters
coef(fit)

# Get confidence intervals
head(confint(fit))

# For regression: get fitted values
head(fitted(fit_reg))

