# rbbnp 0.2.0

## Major Improvements
* Improved performance of `get_est_Ar()` function by approximately 50%, significantly accelerating main function execution time
* Added unit tests to ensure package reliability and stability

## Performance Optimization
* Optimized `get_est_Ar()` function by:
  - Precomputing `avg_phi_log` values for all frequencies
  - Enhancing the optimization algorithm

## Bug Fixes and Edge Case Handling
* Modified `get_sigma()` and `get_sigma_yx()` to handle negative density estimates by returning 0 instead of NaN
* Enhanced `get_conditional_var()` to ensure non-negative variance estimates
* Fixed `sinc()` function to properly handle x = 0 cases

## Function Refactoring
* Added `create_biasBound_config()` function for centralized configuration management
* Added `create_kernel_functions()` factory for flexible kernel function creation
