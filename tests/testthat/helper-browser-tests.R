skip_unless_browser_tests_enabled <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("METABORX_RUN_BROWSER_TESTS"), "true"),
    "Set METABORX_RUN_BROWSER_TESTS=true to run browser-based shinytest2 tests."
  )
}
