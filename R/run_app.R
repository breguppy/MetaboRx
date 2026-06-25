#' Launches MetaboRx Shiny App
#' @export

run_app <- function() {
  check_required_dependencies()
  options(shiny.launch.browser = TRUE, shiny.testmode = FALSE)
  shiny::shinyApp(ui = app_ui(), server = app_server)
}
