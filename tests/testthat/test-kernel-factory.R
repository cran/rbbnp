test_that("create_kernel_functions works with different options", {
  # Test different kernel types
  kernels <- c("Schennach2004", "sinc", "normal", "epanechnikov")
  
  for(k in kernels) {
    result <- create_kernel_functions(kernel.fun = k)
    expect_type(result, "list")
    expect_true(all(c("kernel", "kernel_ft") %in% names(result)))
    expect_type(result$kernel, "closure")
    expect_type(result$kernel_ft, "closure")
  }
}) 