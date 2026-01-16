#' Export module
#' @keywords internal
#' @noRd

mod_export_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "4. Export All", 
    value = "tab_export",
    card(
    card_title(tags$h4("Download all to get a ", icon("folder"), " zipped folder containing:")),
    fluidRow(
      column(4, tags$h5(icon("file-excel"), "Excel files:"),
             tags$ul(style = "list-style-type: none;",
               tags$li(icon("file-excel"),"Corrected and Transformed Data"),
               tags$li(icon("file-excel"),"Missing Value Summary "),
               tags$li(icon("file-excel"),"RSD Summaries"),
               tags$ul(style = "list-style-type:none;",
                       tags$li("Corrected Data RSD Summary"),
                       tags$li("Transformed Data RSD Summary")),
               tags$li(icon("file-excel"),"Extreme Value Summary"),
               tags$li(icon("file-excel"),"Metabolite Correlations")
             )),
      column(4, tags$h5(icon("folder"), " figures"),
             tags$ul(style = "list-style-type: none;",
                     tags$li(icon("folder"), " metabolite figures"),
                     tags$li(icon("folder"), " RSD figures"),
                     tags$li(icon("folder"), "PCA plots")
             )),
      column(4, tags$h5(icon("file-circle-check"), " quality_report.html"),
             tags$ul(style = "list-style-type: none;",
                     tags$li("Report describing all summaries, preprocessing steps, and figures generated from the app.")
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
        temp_params <- p()
        temp_params$rsd_compare <- "filtered_cor_data"
        stats_xlsx_path1 <- file.path(base_dir, sprintf("rsd_stats_%s_%s.xlsx", temp_params$rsd_compare, Sys.Date()))
        stats_wb1 <- export_stats_xlsx(temp_params, d())
        openxlsx::saveWorkbook(stats_wb1, stats_xlsx_path1, overwrite = TRUE)
        
        temp_params$rsd_compare <- "transformed_cor_data"
        stats_xlsx_path2 <- file.path(base_dir, sprintf("rsd_stats_%s_%s.xlsx", temp_params$rsd_compare, Sys.Date()))
        stats_wb2 <- export_stats_xlsx(temp_params, d())
        openxlsx::saveWorkbook(stats_wb2, stats_xlsx_path2, overwrite = TRUE)
        
        # Create and save metabolite correlations data file
        corr_xlsx_path <- file.path(base_dir, sprintf("metabolite_correlations_%s.xlsx", Sys.Date()))
        corr_wb <- export_corr_xlsx(d()$all_corr)
        openxlsx::saveWorkbook(corr_wb, corr_xlsx_path, overwrite = TRUE)
        
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
          basename(stats_xlsx_path1),
          basename(outlier_xlsx_path),
          basename(corr_xlsx_path),
          "quality_report.html",
          if (!identical(p()$transform, "none")) basename(stats_xlsx_path2)
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