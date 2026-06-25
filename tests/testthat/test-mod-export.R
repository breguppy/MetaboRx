testthat::skip_on_cran()
skip_unless_browser_tests_enabled()
testthat::skip_if_not_installed("shinytest2")
testthat::skip_if_not_installed("chromote")

test_that("export module creates a zip bundle", {
  options(shiny.port = NULL, shiny.launch.browser = FALSE)
  Sys.unsetenv("SHINY_PORT")

  app <- shinytest2::AppDriver$new(
    test_path("_apps/mod_export"),
    name = "mod_export_basic",
    seed = 123,
    variant = shinytest2::platform_variant(),
    shiny_args = list(host = "127.0.0.1", port = httpuv::randomPort()),
    view = "none",
    load_timeout = 20000
  )
  on.exit(app$stop(), add = TRUE)

  # Wait for the export UI to exist / be rendered
  app$wait_for_value(output = "export-download_all_ui")
  app$wait_for_idle()

  # Trigger and download the zip bundle
  if ("download" %in% ls(app)) {
    zip_path <- app$download("export-download_all_zip")
  } else {
    zip_path <- app$get_download("export-download_all_zip")
  }

  expect_true(file.exists(zip_path))
  expect_gt(file.size(zip_path), 0)

  lst <- utils::unzip(zip_path, list = TRUE)
  expect_true(any(grepl("^figures/?", lst$Name)))
  # expect_true(any(grepl("^missing_value_counts_.*\\.xlsx$", lst$Name)))
  expect_true(any(grepl("^corrected_data_.*\\.xlsx$", lst$Name)))
  expect_true(any(grepl("^rsd_stats_.*\\.xlsx$", lst$Name)))
  expect_true(any(grepl("^extreme_values_.*\\.xlsx$", lst$Name)))
  expect_true(any(grepl("^quality_report\\.html$", lst$Name)))
})
