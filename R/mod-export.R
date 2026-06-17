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
        column(
          4,
          tags$h5(icon("file-excel"), " Excel files:"),
          tags$ul(
            style = "list-style-type: none;",
            tags$li(icon("file-excel"), "Corrected and Transformed Data"),
            tags$li(icon("file-excel"), "Missing Value Summary"),
            tags$li(icon("file-excel"), "RSD Summaries"),
            tags$ul(
              style = "list-style-type:none;",
              tags$li("Corrected Data RSD Summary"),
              tags$li("Transformed Data RSD Summary")
            ),
            tags$li(icon("file-excel"), "Extreme Value Summary"),
            tags$li(icon("file-excel"), "Metabolite Correlations (optional)"),
            tags$li(icon("file-excel"), "PCA loadings workbook (inside PCA plots folder)")
          )
        ),
        column(
          4,
          tags$h5(icon("folder"), " Figures"),
          tags$ul(
            style = "list-style-type: none;",
            tags$li(icon("folder"), " Metabolite figures"),
            tags$li(icon("folder"), " RSD figures"),
            tags$li(icon("folder"), " PCA plots")
          )
        ),
        column(
          4,
          tags$h5(icon("file-circle-check"), " quality_report.html"),
          tags$ul(
            style = "list-style-type: none;",
            tags$li("Report describing all summaries, preprocessing steps, and figures generated from the app.")
          )
        )
      ),
      uiOutput(ns("download_all_ui"))
    )
  )
}

.export_bundle_cache <- function() {
  new.env(parent = emptyenv())
}

.export_cached_rsd_plot_data <- function(cache, compare_to, p, d) {
  compare_to <- compare_to %||% "filtered_cor_data"

  if (is.null(cache[["rsd_plot_data"]])) {
    cache[["rsd_plot_data"]] <- list()
  }

  if (is.null(cache[["rsd_plot_data"]][[compare_to]])) {
    temp_params <- p
    temp_params$rsd_compare <- compare_to
    cache[["rsd_plot_data"]][[compare_to]] <- .get_rsd_plot_data(temp_params, d)
  }

  cache[["rsd_plot_data"]][[compare_to]]
}

.export_cached_rsd_results <- function(cache, compare_to, p, d) {
  .export_cached_rsd_plot_data(cache, compare_to, p, d)$rsd_results
}

.export_cached_rsd_plots <- function(cache, p, d) {
  if (is.null(cache[["rsd_plots"]])) {
    cache[["rsd_plots"]] <- make_all_rsd_plots(
      p,
      d,
      rsd_cache = cache[["rsd_plot_data"]] %||% list()
    )
  }

  cache[["rsd_plots"]]
}

.export_cached_pca_plots <- function(cache, p, d) {
  if (is.null(cache[["pca_plots"]])) {
    cache[["pca_plots"]] <- suppressWarnings(make_all_pca_plots(p, d, d$cleaned$meta_df))
  }

  cache[["pca_plots"]]
}

.export_cached_pca_pair <- function(cache, p, d) {
  pca_plots <- .export_cached_pca_plots(cache, p, d)
  pca_compare <- p$pca_compare %||% "filtered_cor_data"
  pca_pair <- pca_plots$pca_pairs[[pca_compare]]

  if (is.null(pca_pair)) {
    return(NULL)
  }

  pca_pair
}

.export_cached_hotelling <- function(cache, p, d) {
  if (is.null(cache[["hotelling_res"]])) {
    cache[["hotelling_res"]] <- d$hotelling_res
    if (is.null(cache[["hotelling_res"]])) {
      cache[["hotelling_res"]] <- detect_hotelling_nonqc_dual_z(d$filtered_corrected$df_no_mv, p)
    }
  }

  cache[["hotelling_res"]]
}

mod_export_server <- function(id, data, params) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    d <- reactive(data())
    p <- reactive(params())
    
    output$download_all_ui <- renderUI({
      req(d()$filtered, d()$filtered_corrected)
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 800px; display: inline-block;",
          downloadButton(
            outputId = ns("download_all_zip"),
            label = "Download All",
            class = "btn btn-secondary btn-lg"
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
        export_date <- Sys.Date()
        export_params <- p()
        export_data <- d()
        export_cache <- .export_bundle_cache()

        base_dir <- tempfile("bundle_")
        dir.create(base_dir)
        on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)

        mv_xlsx_path <- file.path(base_dir, sprintf("missing_value_counts_%s.xlsx", export_date))
        mv_wb <- export_mv_xlsx(export_params, export_data)
        openxlsx::saveWorkbook(mv_wb, mv_xlsx_path, overwrite = TRUE)

        xlsx_path <- file.path(base_dir, sprintf("corrected_data_%s.xlsx", export_date))
        wb <- export_xlsx(export_params, export_data)
        openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

        temp_params <- export_params
        temp_params$rsd_compare <- "filtered_cor_data"
        stats_xlsx_path1 <- file.path(
          base_dir,
          sprintf("rsd_stats_%s_%s.xlsx", temp_params$rsd_compare, export_date)
        )
        stats_wb1 <- export_stats_xlsx(
          temp_params,
          export_data,
          rsd_results = .export_cached_rsd_results(
            export_cache,
            "filtered_cor_data",
            export_params,
            export_data
          )
        )
        openxlsx::saveWorkbook(stats_wb1, stats_xlsx_path1, overwrite = TRUE)

        if (!identical(export_params$transform, "none")) {
          temp_params$rsd_compare <- "transformed_cor_data"
          stats_xlsx_path2 <- file.path(
            base_dir,
            sprintf("rsd_stats_%s_%s.xlsx", temp_params$rsd_compare, export_date)
          )
          stats_wb2 <- export_stats_xlsx(
            temp_params,
            export_data,
            rsd_results = .export_cached_rsd_results(
              export_cache,
              "transformed_cor_data",
              export_params,
              export_data
            )
          )
          openxlsx::saveWorkbook(stats_wb2, stats_xlsx_path2, overwrite = TRUE)
        }

        if (!is.null(export_data$all_corr)) {
          corr_xlsx_path <- file.path(base_dir, sprintf("metabolite_correlations_%s.xlsx", export_date))
          corr_wb <- export_corr_xlsx(export_data$all_corr)
          openxlsx::saveWorkbook(corr_wb, corr_xlsx_path, overwrite = TRUE)
        }

        outlier_xlsx_path <- file.path(base_dir, sprintf("extreme_values_%s.xlsx", export_date))
        outlier_wb <- export_outliers_xlsx(
          export_params,
          export_data,
          hotelling_res = .export_cached_hotelling(export_cache, export_params, export_data)
        )
        openxlsx::saveWorkbook(outlier_wb, outlier_xlsx_path, overwrite = TRUE)

        figs <- export_figures(
          export_params,
          export_data,
          out_dir = base_dir,
          rsd_export = .export_cached_rsd_plots(export_cache, export_params, export_data),
          pca_export = .export_cached_pca_plots(export_cache, export_params, export_data)
        )

        rr <- render_report(
          export_params,
          export_data,
          out_dir = base_dir,
          rsd_plot_data = .export_cached_rsd_plot_data(
            export_cache,
            export_params$rsd_compare,
            export_params,
            export_data
          ),
          pca_pair = .export_cached_pca_pair(export_cache, export_params, export_data),
          hotelling_res = .export_cached_hotelling(export_cache, export_params, export_data)
        )

        rel <- c(
          "figures",
          basename(mv_xlsx_path),
          basename(xlsx_path),
          basename(stats_xlsx_path1),
          basename(outlier_xlsx_path),
          if (!is.null(export_data$all_corr)) basename(corr_xlsx_path),
          "quality_report.html",
          if (!identical(export_params$transform, "none")) basename(stats_xlsx_path2)
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
