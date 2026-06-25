test_that("all user-facing runtime packages are hard dependencies", {
  imports <- .metaborx_imports()

  runtime_packages <- c(
    "randomForest",
    "ggtext",
    "cowplot",
    "zip",
    "pmp"
  )

  expect_setequal(intersect(runtime_packages, imports), runtime_packages)
})

test_that("every runtime namespace is declared as a hard dependency", {
  runtime_packages <- c(
    "shiny",
    "ggplot2",
    "dplyr",
    "tidyr",
    "purrr",
    "tibble",
    "stringr",
    "bslib",
    "shinyjs",
    "shinycssloaders",
    "rlang",
    "htmltools",
    "readxl",
    "openxlsx",
    "rmarkdown",
    "impute",
    "randomForest",
    "ggtext",
    "cowplot",
    "zip",
    "pmp"
  )

  expect_setequal(.metaborx_imports(), runtime_packages)
})

test_that("dependency status reports missing packages and Pandoc", {
  status <- .dependency_status(
    imports = c("shiny", "pmp"),
    namespace_available = function(package) identical(package, "shiny"),
    pandoc_available = function() FALSE
  )

  expect_identical(status$missing_packages, "pmp")
  expect_false(status$pandoc_available)
  expect_false(status$ready)
})

test_that("dependency errors tell beginners how to repair the installation", {
  status <- list(
    missing_packages = c("pmp", "zip"),
    pandoc_available = FALSE,
    ready = FALSE
  )

  message <- .dependency_error_message(status)

  expect_match(message, "MetaboRx is not ready to start", fixed = TRUE)
  expect_match(message, "pmp, zip", fixed = TRUE)
  expect_match(message, "BiocManager::install", fixed = TRUE)
  expect_match(message, "RStudio Desktop", fixed = TRUE)
})

test_that("run_app performs dependency preflight before creating the app", {
  testthat::local_mocked_bindings(
    check_required_dependencies = function() {
      stop("preflight ran", call. = FALSE)
    }
  )

  expect_error(run_app(), "preflight ran", fixed = TRUE)
})
