options(shiny.testmode = TRUE)

library(shiny)
library(shinyjs)
if (!requireNamespace("MetaboRx", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload needed")
  pkgload::load_all(path = "../../..", helpers = FALSE, quiet = TRUE)
}
library(MetaboRx)

# tiny in-memory dataset (no file IO)
df <- data.frame(
  sample = c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12"),
  batch  = c("B1","B1","B1","B1","B1","B1","B1","B1","B1","B1","B1","B1"),
  class  = c("QC","Sample","QC","Sample","QC","Sample","QC","Sample","QC","Sample","QC", "QC"),
  order  = c(1,2,3,4,5,6,7,8,9,10,11,12),
  A = c(10, NA, 12, 13, 11, 14, 16, 10, 9, 15, 12, 11),
  B = c(20, 25, NA, 20, 21, 24, 26, 24, 25, 21, 23, 22),
  check.names = FALSE
)

data_stub <- reactiveVal(list(
  cleaned  = list(df = df),
  filtered = list(df = df)
))
params_stub <- reactiveVal(list(
  sample_col = "sample",
  batch_col  = "batch",
  class_col  = "class",
  order_col  = "order",
  Frule      = NA
))

ui <- fluidPage(
  useShinyjs(),
  tabsetPanel(id = "main_steps",
              MetaboRx:::mod_correct_ui("correct"),
              tabPanel(title = "3. Visualize", value = "tab_visualize", "ok")
  ),
  # invisible test outputs
  verbatimTextOutput("correct_params_json"),
  verbatimTextOutput("transformed_n")
)

server <- function(input, output, session) {
  mod <- MetaboRx:::mod_correct_server(
    "correct",
    data   = function() data_stub(),
    params = function() params_stub()
  )
  
  output$correct_params_json <- renderText({
    p <- mod$params()
    if (is.null(p)) return("null")
    jsonlite::toJSON(p, auto_unbox = TRUE, null = "null")
  })
  output$transformed_n <- renderText({
    tr <- mod$transformed()
    if (is.null(tr)) return(NA_integer_)
    df <- tr$df_no_mv
    if (is.null(df)) df <- tr$df_mv
    if (is.null(df)) return(NA_integer_)
    as.character(nrow(df))
  })
  outputOptions(output, "correct_params_json", suspendWhenHidden = FALSE)
  outputOptions(output, "transformed_n",      suspendWhenHidden = FALSE)
}

shinyApp(ui, server)
