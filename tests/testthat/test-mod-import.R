library(testthat)
library(shinytest2)

testthat::skip_on_cran()
testthat::skip_if_not_installed("shinytest2")
testthat::skip_if_not_installed("chromote")

test_that("import module loads, selects columns, filters, and navigates", {
  # ensure no fixed port is forced
  options(shiny.port = NULL)
  Sys.unsetenv("SHINY_PORT")

  app <- AppDriver$new(
    test_path("_apps/mod_import"),
    name = "mod_import_basic",
    seed = 123,
    variant = platform_variant(),
    options = list(shiny.testmode = TRUE),
    shiny_args = list(host = "127.0.0.1", port = httpuv::randomPort()), # not 0
    view = "none",
    load_timeout = 20000
  )

  options(shiny.port = NULL)
  Sys.unsetenv("SHINY_PORT")

  options(shiny.launch.browser = FALSE)

  csv <- test_path("fixtures/raw_small.csv")
  app$upload_file("import-file1" = csv)

  app$wait_for_value(output = "import-contents")
  app$wait_for_js("document.getElementById('import-sample_col') !== null")

  app$set_inputs(
    "import-sample_col" = "sample",
    "import-batch_col"  = "batch",
    "import-class_col"  = "class",
    "import-order_col"  = "order"
  )
  app$wait_for_js("Shiny.shinyapp.$inputValues['import-order_col'] === 'order'")

  app$set_inputs("import-withhold_cols" = TRUE)
  app$set_inputs("import-n_withhold" = 1)
  app$wait_for_js("document.getElementById('import-withhold_col_1') !== null")
  app$set_inputs("import-withhold_col_1" = "metabolite_A")

  app$wait_for_js(
    "document.getElementById('import-mv_cutoff') !== null",
    timeout = 10000
  )
  app$set_inputs("import-mv_cutoff" = 20, wait_ = FALSE)
  app$wait_for_value(output = "import-basic_info")
  app$wait_for_value(output = "import-filter_info")
  app$wait_for_idle()

  # Read outputs instead of exports
  app$wait_for_value(output = "params_json")
  params <- jsonlite::fromJSON(app$get_value(output = "params_json"))
  cleaned_n <- as.integer(app$get_value(output = "cleaned_n"))

  expect_type(params, "list")
  expect_equal(params$sample_col, "sample")
  expect_true("metabolite_A" %in% params$withheld_cols)
  expect_gt(cleaned_n, 0)

  app$wait_for_js(
    "document.getElementById('import-next_correction') !== null",
    timeout = 10000
  )
  app$click("import-next_correction")

  # Wait for Shiny to update the tabset input and assert
  app$wait_for_value(input = "main_steps")
  expect_equal(app$get_value(input = "main_steps"), "tab_correct")
  app$wait_for_js(
    "Shiny.shinyapp && Shiny.shinyapp.$inputValues['main_steps'] === 'tab_correct'",
    timeout = 10000
  )
})
