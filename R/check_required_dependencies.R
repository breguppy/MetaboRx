.metaborx_imports <- function(description = utils::packageDescription("MetaboRx")) {
  imports <- description[["Imports"]]
  if (is.null(imports) || is.na(imports) || !nzchar(imports)) {
    return(character(0))
  }

  imports |>
    strsplit(",") |>
    unlist(use.names = FALSE) |>
    trimws() |>
    sub(pattern = "\\s*\\(.*", replacement = "") |>
    unique()
}

.dependency_status <- function(
  imports,
  namespace_available = function(package) {
    requireNamespace(package, quietly = TRUE)
  },
  pandoc_available = function() {
    rmarkdown::pandoc_available()
  }
) {
  installed <- vapply(imports, namespace_available, logical(1))
  missing_packages <- imports[!installed]
  has_rmarkdown <- !("rmarkdown" %in% imports) || isTRUE(installed[["rmarkdown"]])
  has_pandoc <- has_rmarkdown && isTRUE(pandoc_available())

  list(
    missing_packages = missing_packages,
    pandoc_available = has_pandoc,
    ready = length(missing_packages) == 0L && has_pandoc
  )
}

.dependency_error_message <- function(status) {
  problems <- "MetaboRx is not ready to start."

  if (length(status$missing_packages) > 0L) {
    problems <- c(
      problems,
      paste0(
        "Missing R packages: ",
        paste(status$missing_packages, collapse = ", "),
        "."
      ),
      "Repair the installation by copying these commands into the RStudio Console:",
      "install.packages(c('BiocManager', 'remotes'))",
      "BiocManager::install(c('impute', 'pmp'), ask = FALSE, update = FALSE)",
      paste0(
        "remotes::install_github('breguppy/MetaboRx', ",
        "dependencies = NA, upgrade = 'never')"
      )
    )
  }

  if (!isTRUE(status$pandoc_available)) {
    problems <- c(
      problems,
      paste(
        "Pandoc was not found. Install or update RStudio Desktop, restart",
        "RStudio, and try again. RStudio Desktop includes Pandoc, which",
        "MetaboRx uses to create the HTML quality report."
      )
    )
  }

  paste(problems, collapse = "\n")
}

#' Check required MetaboRx dependencies
#'
#' Checks that MetaboRx is ready to run
#'
#' @details
#' This function derives the runtime package list from the package's
#' \code{Imports} field, checks that every package is installed, and verifies
#' that Pandoc is available for the HTML quality report included in Download
#' All. Core dependencies are installed automatically with MetaboRx. RStudio
#' Desktop includes Pandoc.
#'
#' @return
#' Invisibly returns a list containing missing packages, Pandoc availability,
#' and an overall readiness flag. Stops with repair instructions when the
#' installation is incomplete.
#'
#' @examples
#' \dontrun{
#' MetaboRx::check_required_dependencies()
#' }
#' @export
check_required_dependencies <- function() {
  status <- .dependency_status(.metaborx_imports())

  if (!status$ready) {
    stop(.dependency_error_message(status), call. = FALSE)
  }

  invisible(status)
}
