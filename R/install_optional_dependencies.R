#' Install optional MetaboRx dependencies
#'
#' Installs optional packages used to develop and test MetaboRx.
#'
#' @details
#' Normal MetaboRx users do not need to call this function. Every package used
#' by the app's analysis and export features is listed in \code{Imports} and is
#' installed with MetaboRx. This helper installs packages used by contributors
#' for automated tests, browser tests, benchmarks, and package development.
#'
#' @return
#' Invisibly returns a list with character vectors of missing CRAN and
#' Bioconductor packages that were installed.
#'
#' @examples
#' \dontrun{
#' MetaboRx::install_optional_dependencies()
#' }
#' @export
install_optional_dependencies <- function() {
  cran_pkgs <- c(
    "httpuv",
    "jsonlite",
    "testthat",
    "shinytest2",
    "chromote",
    "pkgload"
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
    message("All MetaboRx development dependencies are already installed.")
  }

  invisible(list(
    cran = cran_missing,
    bioc = character(0)
  ))
}
