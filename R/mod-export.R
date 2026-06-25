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
    d <- reactive(data())
    p <- reactive(params())

    output$download_all_ui <- renderUI({
      req(d()$filtered, d()$filtered_corrected)
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 800px; display: inline-block;",
          downloadButton(
            outputId = session$ns("download_all_zip"),
            label = "Download All",
            class = "btn btn-secondary btn-lg"
          )
        )
      )
    })

    output$download_all_zip <- downloadHandler(
      filename = function() {
        sprintf("corrected_data_plots_report_%s.zip", Sys.Date())
      },
      content = function(file) {
        export_bundle(p = p(), d = d(), file = file)
      }
    )
  })
}
