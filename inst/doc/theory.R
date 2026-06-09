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

## ----fourier-demo-------------------------------------------------------------
# Generate sample data
X <- gen_sample_data(size = 500, dgp = "2_fold_uniform", seed = 42)

# Estimate density
fit <- biasBound_density(X, h = 0.08, kernel.fun = "Schennach2004")

# View detected smoothness parameters
coef(fit)

## ----ft-plot------------------------------------------------------------------
# Visualize Fourier transform fit
plot(fit, type = "ft")

## ----ci-visualization---------------------------------------------------------
# The plot shows both bands
plot(fit)

## ----kernel-comparison, fig.width=6, fig.height=10.5, out.width="100%"--------
library(gridExtra)

fit_sch <- biasBound_density(X, kernel.fun = "Schennach2004")
fit_sinc <- biasBound_density(X, kernel.fun = "sinc")

grid.arrange(
  plot(fit_sch) + ggtitle("Schennach2004 (recommended)"),
  plot(fit_sinc) + ggtitle("Sinc kernel"),
  ncol = 1
)

## ----regression-theory--------------------------------------------------------
# Generate regression data
Y <- sin(2 * pi * X) + rnorm(500, sd = 0.3)

# Estimate with bias bounds
fit_reg <- biasBound_condExpectation(Y, X, h = 0.1)

# View smoothness parameters
coef(fit_reg)

## ----bw-selection-------------------------------------------------------------
h_cv <- select_bandwidth(X, method = "cv", kernel.fun = "Schennach2004")
h_silv <- select_bandwidth(X, method = "silverman", kernel.fun = "normal")

cat("CV bandwidth:", round(h_cv, 4), "\n")
cat("Silverman bandwidth:", round(h_silv, 4))

## ----optimal-vs-under, fig.width=6, fig.height=10.5, out.width="100%"---------
result_opt <- biasBound_density(X, h = h_cv, kernel.fun = "Schennach2004")
result_under <- biasBound_density(X, h = h_cv * 0.5, kernel.fun = "Schennach2004")

grid.arrange(
  plot(result_opt) + ggtitle(paste0("Optimal bandwidth (h = ", round(h_cv, 3), ")")),
  plot(result_under) + ggtitle(paste0("Undersmoothed (h = ", round(h_cv/2, 3), ")")),
  ncol = 1
)

