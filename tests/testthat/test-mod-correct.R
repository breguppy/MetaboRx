testthat::skip_on_cran()
skip_unless_browser_tests_enabled()
testthat::skip_if_not_installed("shinytest2")
testthat::skip_if_not_installed("chromote")

test_that("correct module runs and navigates", {
  options(shiny.port = NULL, shiny.launch.browser = FALSE)
  Sys.unsetenv("SHINY_PORT")

  app <- shinytest2::AppDriver$new(
    test_path("_apps/mod_correct"),
    name = "mod_correct_basic",
    seed = 123,
    variant = shinytest2::platform_variant(),
    shiny_args = list(host = "127.0.0.1", port = httpuv::randomPort()),
    view = "none",
    load_timeout = 20000
  )
  on.exit(app$stop(), add = TRUE)

  # 1) Wait for the rendered UI chunks, not the inner inputs
  app$wait_for_value(output = "correct-qcImpute")
  app$wait_for_value(output = "correct-sampleImpute")
  app$wait_for_value(output = "correct-correctionMethod")

  # 2) Resolve the actual select IDs from those containers
  qc_id <- app$get_js("let s=document.querySelector('#correct-qcImpute select'); s ? s.id : null")
  sam_id <- app$get_js("let s=document.querySelector('#correct-sampleImpute select'); s ? s.id : null")
  cor_id <- app$get_js("let s=document.querySelector('#correct-correctionMethod select'); s ? s.id : null")

  # 3) Pick a valid value for each select (first non-empty option)
  qc_val <- app$get_js("let s=document.querySelector('#correct-qcImpute select'); if(!s) null; else { for (const o of s.options){ if(o.value && !o.disabled) return o.value } return s.options.length? s.options[0].value : null }")
  sam_val <- app$get_js("let s=document.querySelector('#correct-sampleImpute select'); if(!s) null; else { for (const o of s.options){ if(o.value && !o.disabled) return o.value } return s.options.length? s.options[0].value : null }")
  cor_val <- app$get_js("let s=document.querySelector('#correct-correctionMethod select'); if(!s) null; else { for (const o of s.options){ if(o.value && !o.disabled) return o.value } return s.options.length? s.options[0].value : null }")

  # 4) Set inputs using the resolved IDs
  if (!is.null(qc_id) && !is.null(qc_val)) do.call(app$set_inputs, setNames(list(qc_val), qc_id))
  if (!is.null(sam_id) && !is.null(sam_val)) do.call(app$set_inputs, setNames(list(sam_val), sam_id))
  if (!is.null(cor_id) && !is.null(cor_val)) do.call(app$set_inputs, setNames(list(cor_val), cor_id))

  # 5) Run correction
  app$click("correct-correct")

  # 6) Wait for downstream values exposed by the harness
  app$wait_for_value(output = "correct-cor_data")
  app$wait_for_value(output = "correct-post_cor_filter_info")
  app$wait_for_value(output = "correct_params_json")
  app$wait_for_value(output = "transformed_n")

  # 7) Assertions
  params <- jsonlite::fromJSON(app$get_value(output = "correct_params_json"))
  n_tr <- as.integer(app$get_value(output = "transformed_n"))
  expect_type(params, "list")
  expect_gt(n_tr, 0)

  # 8) Navigate to visualize tab
  app$click("correct-next_visualization")
  app$wait_for_value(input = "main_steps")
  expect_equal(app$get_value(input = "main_steps"), "tab_visualize")
})
