#' Create kernel functions based on configuration
#'
#' @param kernel.fun A string specifying the kernel function to be used.
#' @param if_approx_kernel Logical. If TRUE, uses approximations for the kernel function.
#' @param kernel.resol The resolution for kernel function approximation.
#' @param X Optional data vector used to determine approximation range.
#' @param h Optional bandwidth parameter used for approximation range.
#'
#' @return A list containing kernel function and its Fourier transform
#'
#' @export
create_kernel_functions <- function(kernel.fun = "Schennach2004", 
                                   if_approx_kernel = TRUE, 
                                   kernel.resol = 1000,
                                   X = NULL,
                                   h = 0.09) {
  
  if (kernel.fun == "Schennach2004") {
    if (if_approx_kernel) {
      if (!is.null(X) && !is.null(h)) {
        # determine the range of the interpolation
        u_lb <- -(max(X) - min(X)) / h * 10
        u_ub <- (max(X) - min(X)) / h * 10
        
        inf_k <- function(u) fun_approx(u, u_lb = u_lb, u_ub = u_ub, resol = kernel.resol)
      } else {
        # Default range if X and h not provided
        inf_k <- function(u) fun_approx(u, u_lb = -100, u_ub = 100, resol = kernel.resol)
      }
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
  
  return(list(kernel = inf_k, kernel_ft = inf_k_ft))
} 