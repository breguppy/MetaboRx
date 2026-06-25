testthat::skip_on_cran()
skip_unless_browser_tests_enabled()
testthat::skip_if_not_installed("shinytest2")
testthat::skip_if_not_installed("chromote")

test_that("visualize module renders and navigates", {
  options(shiny.port = NULL, shiny.launch.browser = FALSE)
  Sys.unsetenv("SHINY_PORT")

  app <- shinytest2::AppDriver$new(
    test_path("_apps/mod_visualize"),
    name = "mod_visualize_basic",
    seed = 123,
    variant = shinytest2::platform_variant(),
    shiny_args = list(host = "127.0.0.1", port = httpuv::randomPort()),
    view = "none",
    load_timeout = 20000
  )
  on.exit(app$stop(), add = TRUE)

  # wait for metabolite selector to render, then choose first valid option
  app$wait_for_js("document.querySelector('#visualize-met_plot_selectors select') !== null")
  met_id <- app$get_js("let s=document.querySelector('#visualize-met_plot_selectors select'); s ? s.id : null")
  met_val <- app$get_js("let s=document.querySelector('#visualize-met_plot_selectors select'); if(!s) null; for (const o of s.options){ if(o.value && !o.disabled) return o.value } null")
  if (!is.null(met_id) && !is.null(met_val)) do.call(app$set_inputs, setNames(list(met_val), met_id))

  # discover other selects created by ui_rsd_eval/ui_pca_eval and set first valid choices
  ids <- app$get_js("Array.from(document.querySelectorAll('select[id^=\"visualize-\"]')).map(e=>e.id)")
  for (sid in ids) {
    val <- app$get_js(sprintf(
      "let s=document.getElementById('%s'); if(!s) null; for (const o of s.options){ if(o.value && !o.disabled) return o.value } null",
      sid
    ))
    if (!is.null(val)) do.call(app$set_inputs, setNames(list(val), sid))
  }

  # plots should render
  app$wait_for_value(output = "visualize-metab_scatter")
  app$wait_for_value(output = "visualize-rsd_comparison_plot")
  app$wait_for_value(output = "visualize-pca_plot")
  app$wait_for_value(output = "visualize-pca_loading_plot")

  # navigate to export tab
  app$click("visualize-next_export")
  app$wait_for_value(input = "main_steps")
  expect_equal(app$get_value(input = "main_steps"), "tab_export")
})
