.onLoad <- function(libname, pkgname) {
  # Define the PATH variable
  DATA_PATH <- system.file("data", package = pkgname)
  # Assign the DATA_PATH variable to the package's namespace
  assign("DATA_PATH", DATA_PATH, envir = parent.env(environment()))

  # Define the path variable for external (non-R) data files
  EXT_DATA_PATH <- system.file("extdata", package = pkgname)
  # Assign the EXT_DATA_PATH variable to the package's namespace
  assign("EXT_DATA_PATH", EXT_DATA_PATH, envir = parent.env(environment()))
}

#' The Path to the Data Folder
#'
#' This variable provides the path to the `data` folder within the package.
#' @name DATA_PATH
#' @return The path to the package's internal data folder as a character string.
#' @export
NULL

#' The Path to the External Data Folder for Non-R Data Files
#'
#' This variable provides the path to the `extdata` folder within the package,
#' where non-standard R data files are stored.
#' @name EXT_DATA_PATH
#' @return The path to the package's external data folder (for non-standard R data files) as a character string.
#' @export
NULL
