options(shiny.testmode = TRUE)
library(shiny)
if (!requireNamespace("MetaboRx", quietly = TRUE)) {
  pkgload::load_all(path = "../../..", helpers = FALSE, quiet = TRUE)
}
library(MetaboRx)

shinyApp(
  ui     = MetaboRx:::app_ui(),
  server = function(input, output, session) MetaboRx:::app_server(input, output, session)
)