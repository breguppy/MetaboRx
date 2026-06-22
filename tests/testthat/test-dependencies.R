test_that("all user-facing runtime packages are hard dependencies", {
  description_path <- testthat::test_path("..", "..", "DESCRIPTION")
  description <- read.dcf(description_path)
  imports <- description[1, "Imports"] |>
    strsplit(",") |>
    unlist() |>
    trimws() |>
    sub(pattern = "\\s*\\(.*", replacement = "")

  runtime_packages <- c(
    "randomForest",
    "ggtext",
    "cowplot",
    "zip",
    "pmp",
    "knitr"
  )

  expect_setequal(intersect(runtime_packages, imports), runtime_packages)
})

test_that("every runtime namespace is declared as a hard dependency", {
  package_root <- testthat::test_path("..", "..")
  description <- read.dcf(file.path(package_root, "DESCRIPTION"))
  imports <- description[1, "Imports"] |>
    strsplit(",") |>
    unlist() |>
    trimws() |>
    sub(pattern = "\\s*\\(.*", replacement = "")

  runtime_files <- c(
    list.files(
      file.path(package_root, "R"),
      pattern = "\\.[Rr]$",
      full.names = TRUE,
      recursive = TRUE
    ),
    list.files(
      file.path(package_root, "inst"),
      pattern = "\\.(R|Rmd)$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
  runtime_text <- paste(
    vapply(runtime_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1)),
    collapse = "\n"
  )
  namespace_matches <- regmatches(
    runtime_text,
    gregexpr("[A-Za-z][A-Za-z0-9.]*:::{1,2}", runtime_text, perl = TRUE)
  )[[1]]
  namespace_packages <- unique(sub(":::{1,2}$", "", namespace_matches))
  standard_packages <- c("base", "grDevices", "grid", "stats", "tools", "utils")
  internal_packages <- "MetaboRx"

  undeclared <- setdiff(
    namespace_packages,
    c(imports, standard_packages, internal_packages)
  )
  expect_length(undeclared, 0L)
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
