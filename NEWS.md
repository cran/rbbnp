# rbbnp 0.3.0

## Major Improvements
* Added automatic bandwidth selection using cross-validation or Silverman's rule of thumb
* Improved performance of `get_est_Ar()` function by approximately 50%, significantly accelerating main function execution time
* Improved code structure with centralized configuration management
* Added unit tests

## New Features
* Implemented automatic bandwidth selection with two methods:
  - Cross-validation (`h_method = "cv"`) for optimal accuracy
  - Silverman's rule of thumb (`h_method = "silverman"`) for faster computation
* Added `create_biasBound_config()` function for centralized configuration management
* Added `create_kernel_functions()` factory for flexible kernel function creation
* Changed default bandwidth parameter `h` to `NULL` to enable automatic selection

## Performance Optimization
* Optimized `get_est_Ar()` function by:
  - Precomputing `avg_phi_log` values for all frequencies
  - Enhancing the optimization algorithm for determining smoothness parameters
* Increased default resolution for kernel function approximation

## Bug Fixes and Edge Case Handling
* Modified `get_sigma()` and `get_sigma_yx()` to handle negative density estimates by returning 0 instead of NaN
* Enhanced `get_conditional_var()` to ensure non-negative variance estimates
* Fixed `sinc()` function to properly handle x = 0 cases
* Improved handling of NA values in Fourier transform calculations

## Documentation
* Updated examples to demonstrate automatic bandwidth selection
* Improved parameter descriptions throughout the package
