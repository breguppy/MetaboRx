test_that("full app loads and basic flow works", {
  options(shiny.port = NULL, shiny.launch.browser = FALSE)
  Sys.unsetenv("SHINY_PORT")
  smoke_test_packages <- c("httpuv", "shinytest2")
  missing_packages <- smoke_test_packages[
    !purrr::map_lgl(smoke_test_packages, requireNamespace, quietly = TRUE)
  ]
  if (length(missing_packages) > 0 && identical(Sys.getenv("CI"), "true")) {
    testthat::fail(
      paste(
        "Full app smoke test dependencies are missing in CI:",
        paste(missing_packages, collapse = ", ")
      )
    )
  }
  testthat::skip_if(
    length(missing_packages) > 0,
    paste(
      "Full app smoke test dependencies are not installed:",
      paste(missing_packages, collapse = ", ")
    )
  )

  app <- shinytest2::AppDriver$new(
    test_path("_apps/full_app"),
    name = "full_app_smoke",
    variant = shinytest2::platform_variant(),
    shiny_args = list(host = "127.0.0.1", port = httpuv::randomPort()),
    view = "none",
    load_timeout = 20000
  )
  on.exit(app$stop(), add = TRUE)

  # upload small fixture and drive across tabs
  csv <- test_path("fixtures/raw_small.csv")
  app$upload_file("import-file1" = csv)
  app$wait_for_value(output = "import-contents")

  # select meta columns
  app$wait_for_js("document.getElementById('import-sample_col') !== null")
  app$set_inputs(
    "import-sample_col" = "sample", "import-batch_col" = "batch",
    "import-class_col" = "class", "import-order_col" = "order"
  )
  app$wait_for_value(output = "import-filter_info")
  app$wait_for_value(output = "import-next_correction_ui")
  app$wait_for_js("
    (() => {
      const btn = document.getElementById('import-next_correction');
      return btn !== null && !btn.disabled;
    })()
  ")
  app$click("import-next_correction")
  app$wait_for_js("Shiny.shinyapp && Shiny.shinyapp.$inputValues['main_steps'] === 'tab_correct'")
  expect_equal(app$get_value(input = "main_steps"), "tab_correct")

  # minimal correction
  app$wait_for_value(output = "correct-correctionMethod")
  cor_id <- app$get_js("let s=document.querySelector('#correct-correctionMethod select'); s && s.id")
  cor_val <- app$get_js("let s=document.querySelector('#correct-correctionMethod select'); if(!s) null; (Array.from(s.options).find(o=>o.value))?.value")
  if (!is.null(cor_id) && !is.null(cor_val)) {
    do.call(app$set_inputs, stats::setNames(list(cor_val), cor_id))
  }
  app$click("correct-correct", wait_ = FALSE)
  app$wait_for_value(output = "correct-cor_data")
  app$wait_for_value(output = "correct-correlation_slider")
  app$wait_for_js("
    (() => {
      const btn = document.getElementById('correct-next_visualization');
      return btn !== null && !btn.disabled;
    })()
  ")
  app$click("correct-next_visualization")

  app$wait_for_js("Shiny.shinyapp && Shiny.shinyapp.$inputValues['main_steps'] === 'tab_visualize'")
  expect_equal(app$get_value(input = "main_steps"), "tab_visualize")

  # navigate to export
  app$click("viz-next_export")
  app$wait_for_js("Shiny.shinyapp && Shiny.shinyapp.$inputValues['main_steps'] === 'tab_export'")
  expect_equal(app$get_value(input = "main_steps"), "tab_export")
})
