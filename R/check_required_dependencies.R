#' Check required QCcorrection dependencies
#'
#' Checks that packages required by QCcorrection are installed
#'
#' @details
#' This function checks that core package dependencies listed in \code{Imports}
#' are installed. Core dependencies are installed automatically when
#' QCcorrection itself is installed.
#'
#' @return
#' Invisibly returns a list of missing require packages
#'
#' @examples
#' \dontrun{
#' QCcorrection::check_required_dependencies()
#' }
#' @export
check_required_dependencies <- function() {
  required_pkgs <- c(
    "shiny", "ggplot2", "dplyr", "tidyr", "purrr", "tibble",
    "stringr", "bslib", "shinyjs", "shinycssloaders", "rlang",
    "htmltools", "readxl", "openxlsx", "rmarkdown",
    "impute"
  )

  missing <- required_pkgs[!vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]

  if (length(missing) > 0L) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
