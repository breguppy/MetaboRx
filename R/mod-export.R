#' Export module
#' @keywords internal
#' @noRd

mod_export_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "4. Export All", 
    value = "tab_export",
    card(
    card_title("Download Data, Summaries, Plots, and Report"),
    tags$h4("Download all to get a ", icon("folder"), " zipped folder containing:"),
    fluidRow(
      column(4, tags$h5(icon("file-excel"), "corrected_data_*today's_date*.xlsx"),
             tags$ul(style = "list-style-type: none;",
               tags$li(icon("table"), "0. Raw Data"),
               tags$li(icon("table"), "1. Correction Settings"),
               tags$li(icon("table"), "2. Drift Normalized"),
               tags$li(icon("table"), "3. Scaled or Normalized"),
               tags$li(icon("table"), "4. Grouped Data Organized"),
               tags$li(icon("table"), "5. Grouped Data Fold Change (only when provided a control class)"),
               tags$li(icon("table"), "Appendix1. Metaboanalyst Ready")
             )),
      column(4, tags$h5(icon("file-excel"), "rsd_stats_*today's_date*.xlsx"),
             tags$ul(style = "list-style-type: none;",
               tags$li(icon("table"), "Raw Metabolite RSD"),
               tags$li(icon("table"), "Corrected Metabolite RSD or Transformed Corrected Metabolite RSD"),
               tags$li(icon("table"), "Metabolite RSD Comparison"),
               tags$li(icon("table"), "Raw Class RSD"),
               tags$li(icon("table"), "Corrected Class RSD or Transformed Corrected Class RSD"),
               tags$li(icon("table"), "Class RSD Comparison"),
               )),
      column(4, tags$h5(icon("file-excel"), "missing_value_counts_*today's_date*.xlsx"),
             tags$ul(style = "list-style-type: none;",
               tags$li(icon("table"), "Metabolite"),
               tags$li(icon("table"), "Sample"),
               tags$li(icon("table"), "Class"),
               tags$li(icon("table"), "Batch")
             ))
      ),
      fluidRow(
        column(4, tags$h5(icon("file-excel"), "extreme_values_*today's_date*.xlsx"),
             tags$ul(style = "list-style-type: none;",
               tags$li(icon("table"), "Sample MD"),
               tags$li(icon("table"), "Candidates"),
               tags$li(icon("table"), "Confirmations")
               )),
      column(4, tags$h5(icon("folder"), " figures"),
             tags$ul(style = "list-style-type: none;",
               tags$li(icon("folder"), " metabolite figures"),
               tags$li(icon("folder"), " RSD figures"),
               tags$li(icon("folder"), "PCA plots")
             )),
      column(4, tags$h5(icon("file-pdf"), " correction_report(.pdf/.html)"),
             tags$ul(style = "list-style-type: none;",
               tags$li("Report describing the correction steps and figures for evaluating the correction process.")
             ))
    ),
    uiOutput(ns("download_all_ui"))
  ))
}

mod_export_server <- function(id, data, params) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    d <- reactive(data())
    p <- reactive(params())
    
    output$download_all_ui <- renderUI({
      req(d()$transformed)
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 800px; display: inline-block;",
          downloadButton(
            outputId = ns("download_all_zip"),
            label = "Download All",
            class    = "btn btn-secondary btn-lg"
          )
        ) 
      )
    })
    
    #-- Allow user to download corrected data, figures, and correction report.
    output$download_all_zip <- downloadHandler(
      filename = function() {
        sprintf("corrected_data_plots_report_%s.zip", Sys.Date())
      },
      content = function(file) {
        .require_pkg("zip", "create a zip archive")
        base_dir <- tempfile("bundle_")
        dir.create(base_dir)
        on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)
        
        # create and save missing value count excel
        mv_xlsx_path <- file.path(base_dir, sprintf("missing_value_counts_%s.xlsx", Sys.Date()))
        mv_wb <- export_mv_xlsx(p(), d())
        openxlsx::saveWorkbook(mv_wb, mv_xlsx_path, overwrite = TRUE)
        
        # create and save corrected data file
        xlsx_path <- file.path(base_dir, sprintf("corrected_data_%s.xlsx", Sys.Date()))
        wb <- export_xlsx(p(), d())
        openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
        
        # Create and save rsd stats data file
        stats_xlsx_path <- file.path(base_dir, sprintf("rsd_stats_%s.xlsx", Sys.Date()))
        stats_wb <- export_stats_xlsx(p(), d())
        openxlsx::saveWorkbook(stats_wb, stats_xlsx_path, overwrite = TRUE)
        
        # Create and save outlier data file
        outlier_xlsx_path <- file.path(base_dir, sprintf("extreme_values_%s.xlsx", Sys.Date()))
        outlier_wb <- export_outliers_xlsx(p(), d())
        openxlsx::saveWorkbook(outlier_wb, outlier_xlsx_path, overwrite = TRUE)
        
        # create and save figure folder
        figs <- export_figures(p(), d(), out_dir = base_dir)
        
        # make pdf report
        rr <- render_report(p(), d(), out_dir = base_dir)
        
        # make zip file
        rel <- c(
          "figures",
          basename(mv_xlsx_path),
          basename(xlsx_path),
          basename(stats_xlsx_path),
          basename(outlier_xlsx_path),
          "quality_report.html"
        )
        rel <- rel[file.exists(file.path(base_dir, rel))]
        
        tmpzip <- tempfile(fileext = ".zip")
        on.exit(unlink(tmpzip, force = TRUE), add = TRUE)
        zip::zipr(zipfile = tmpzip, files = rel, root = base_dir)
        
        file.copy(tmpzip, file, overwrite = TRUE)
      }
    )
    
  })
}