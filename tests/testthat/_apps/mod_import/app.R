options(shiny.testmode = TRUE)
library(shiny)

if (!requireNamespace("MetaboRx", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload needed")
  pkgload::load_all(path = "../../..", helpers = FALSE, quiet = TRUE)
}
library(MetaboRx)

ui <- fluidPage(
  tabsetPanel(id = "main_steps",
              MetaboRx:::mod_import_ui("import"),
              tabPanel(title = "2. Correct", value = "tab_correct", "placeholder")  # add this
  )
)

server <- function(input, output, session) {
  mod <- MetaboRx:::mod_import_server("import")
  
  output$params_json <- renderText({
    p <- mod$params()
    if (is.null(p)) return("null")
    jsonlite::toJSON(p, auto_unbox = TRUE, null = "null")
  })
  output$cleaned_n <- renderText({
    cd <- mod$cleaned()
    if (is.null(cd) || is.null(cd$df)) return(NA_integer_)
    as.character(nrow(cd$df))
  })
  
  # force computation even though not in UI
  outputOptions(output, "params_json", suspendWhenHidden = FALSE)
  outputOptions(output, "cleaned_n",   suspendWhenHidden = FALSE)
}

shinyApp(ui, server)