options(shiny.testmode = TRUE)

library(shiny)
if (!requireNamespace("MetaboRx", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload needed")
  pkgload::load_all(path = "../../..", helpers = FALSE, quiet = TRUE)
}
library(MetaboRx)

# minimal synthetic dataset
df <- data.frame(
  sample = c("S1", "S2", "S3", "S4", "S5", "S6"),
  batch = c("B1", "B1", "B1", "B2", "B2", "B2"),
  class = c("QC", "Sample", "Sample", "Sample", "Sample", "QC"),
  order = 1:6,
  A = c(10, 11, 13, 12, 14, 11),
  B = c(20, 22, 21, 23, 24, 21),
  C = c(30, 31, 29, 32, 33, 30),
  check.names = FALSE
)

rv_list <- list(df = df)
data_stub <- reactiveVal(list(
  filtered           = rv_list,
  imputed            = rv_list,
  corrected          = list(df = df, str = "Random Forest"),
  filtered_corrected = list(df_no_mv = df, df_mv = df),
  transformed        = list(df_no_mv = df, df_mv = df),
  cleaned            = list(meta_df = df[, c("sample", "batch", "class", "order")])
))
params_stub <- reactiveVal(list(
  remove_imputed = FALSE,
  transform = "none",
  qcImputeM = "median",
  samImputeM = "median",
  rsd_compare = NULL,
  rsd_cal = NULL,
  pca_compare = NULL,
  color_col = NULL,
  shape_col = NULL,
  fig_format = NULL
))

ui <- fluidPage(
  tabsetPanel(
    id = "main_steps",
    MetaboRx:::mod_visualize_ui("visualize"),
    tabPanel("4. Export", value = "tab_export", "ok")
  )
)

server <- function(input, output, session) {
  MetaboRx:::mod_visualize_server(
    "visualize",
    data   = function() data_stub(),
    params = function() params_stub()
  )
}

shinyApp(ui, server)
