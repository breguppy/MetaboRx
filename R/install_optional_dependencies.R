#' Install optional QCcorrection dependencies
#'
#' Installs optional packages used by QCcorrection for advanced correction,
#' outlier detection, reporting, and testing.
#'
#' @details
#' This function installs packages intended for optional functionality, not the
#' core package dependencies listed in \code{Imports}. Core dependencies are
#' installed automatically when QCcorrection itself is installed.
#'
#' @return
#' Invisibly returns a list with character vectors of missing CRAN and
#' Bioconductor packages that were installed.
#'
#' @examples
#' \dontrun{
#' QCcorrection::install_optional_dependencies()
#' }
#' @export
install_optional_dependencies <- function() {
  cran_pkgs <- c(
    "randomForest",
    "robustbase",
    "outliers",
    "EnvStats",
    "ggtext",
    "httpuv",
    "cowplot",
    "jsonlite",
    "zip",
    "testthat",
    "shinytest2",
    "corpcor",
    "pkgload",
    "knitr"
  )
  
  installed <- rownames(installed.packages())
  cran_missing <- setdiff(cran_pkgs, installed)
  
  if (length(cran_missing) > 0L) {
    install.packages(cran_missing)
    message(
      "Installed CRAN packages: ",
      paste(cran_missing, collapse = ", ")
    )
  }
  
  if (length(cran_missing) == 0L) {
    message("All optional QCcorrection dependencies already installed.")
  }
  
  invisible(list(
    cran = cran_missing,
    bioc = character(0)
  ))
}