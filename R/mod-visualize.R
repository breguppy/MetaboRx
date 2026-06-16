#' @keywords internal


mod_visualize_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "3. Evaluation Metrics and Visualization",
    value = "tab_visualize",
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "3.1 Scatter Plot Evaluation",
        shiny::tags$h6("Visualize correction with metabolite scatter plots"),
        shiny::tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
          shiny::tags$strong("Plot guide"),
          bslib::popover(
            shiny::tags$button(
              type = "button",
              class = "btn btn-link p-0",
              style = "text-decoration:none;",
              shiny::icon("circle-info")
            ),
            report_text_met_scatter(),
            title = "What information does the metabolite scatter plot show?",
            placement = "auto",
            options = list(
              container = "body",
              customClass = "popover-responsive"
            )
          )
        ),
        uiOutput(ns("met_plot_selectors")),
        width = 400
      ),
      plotOutput(ns("metab_scatter"), height = "600px", width = "600px") |> withSpinner(color = "#404040"),
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "3.2 RSD Evaluation",
        tags$h6("Evaluate correction method by the change in relative standard deviation (RSD)."),
        shiny::tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
          shiny::tags$strong("Plot guide"),
          bslib::popover(
            shiny::tags$button(
              type = "button",
              class = "btn btn-link p-0",
              style = "text-decoration:none;",
              shiny::icon("circle-info")
            ),
            report_text_rsd_plots(),
            title = "What information does the RSD comparison plots show?",
            placement = "auto",
            options = list(
              container = "body",
              customClass = "popover-responsive"
            )
          )
        ),
        ui_rsd_eval(ns),
        width = 400
      ),
      plotOutput(ns("rsd_comparison_plot"), height = "540px", width = "900px") |> withSpinner(color = "#404040")
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "3.3 PCA Evaluation",
        tags$h6("Evaluate correction using principal component analysis (PCA)."),
        # info button for PCA plots and loading plots
        shiny::tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
          shiny::tags$strong("Plot guide"),
          bslib::popover(
            shiny::tags$button(
              type = "button",
              class = "btn btn-link p-0",
              style = "text-decoration:none;",
              shiny::icon("circle-info")
            ),
            report_text_pca_plots(),
            title = "What information does the PCA plots and loading plots show?",
            placement = "auto",
            options = list(
              container = "body",
              customClass = "popover-responsive"
            )
          )
        ),
        uiOutput(ns("pca_options")),
        width = 400
      ),
      plotOutput(ns("pca_plot"), height = "530px", width = "1200px") |> withSpinner(color = "#404040"),
      plotOutput(ns("pca_loading_plot"), height = "530px", width = "1000px") |> withSpinner(color = "#404040")
    )),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "3.4 Select Figure Format",
          help = c("Figures will be downloaded in the format selected on the right."),
          width = 400
        ),
        fluidRow(
          column(8, ui_fig_format(ns)),
          column(4, uiOutput(ns("download_fig_zip_btn")))
        ),
        uiOutput(ns("progress_ui"))
      )
    ),
    card(
      actionButton(ns("next_export"), "Next: Export All", class = "btn-primary btn-lg"),
    )
  )
}

mod_visualize_server <- function(id, data, params) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    d <- reactive(data())
    p <- reactive(params())

    #-- Let user select which metabolite to display in scatter plot
    output$met_plot_selectors <- renderUI({
      req(d()$filtered, d()$filtered_corrected)
      df_raw <- d()$filtered$df
      if (isTRUE(p()$remove_imputed)) {
        df_cor <- d()$filtered_corrected$df_mv
      } else {
        df_cor <- d()$filtered_corrected$df_no_mv
      }
      raw_cols <- setdiff(names(df_raw), c("sample", "batch", "class", "order"))
      cor_cols <- setdiff(names(df_cor), c("sample", "batch", "class", "order"))
      cols <- intersect(raw_cols, cor_cols)
      validate(need(length(cols) >= 1, "No overlapping metabolites."))
      selectInput(ns("met_col"), "Metabolite column", choices = cols, selected = cols[1])
    })

    #-- Metabolite scatter plot
    output$metab_scatter <- renderPlot(
      {
        req(input$met_col)
        suppressWarnings({
          print(make_met_scatter(d(), p(), input$met_col))
        })
      },
      res = 120
    )

    #-- RSD comparison plot
    output$rsd_comparison_plot <- renderPlot(execOnResize = FALSE, res = 120, {
      req(input$rsd_compare, input$rsd_cal)

      rsd_data <- rsd_plot_data()
      make_rsd_plot(
        list(
          rsd_compare = input$rsd_compare,
          rsd_cal = input$rsd_cal,
          rsd_plot_type = input$rsd_plot_type,
          remove_imputed = p()$remove_imputed
        ),
        d(),
        rsd_results = rsd_data$rsd_results,
        compared_to = rsd_data$compared_to
      )
    })

    rsd_plot_data <- reactive({
      req(input$rsd_compare)

      .get_rsd_plot_data(
        list(
          rsd_compare = input$rsd_compare,
          remove_imputed = p()$remove_imputed
        ),
        d()
      )
    })

    #-- PCA options
    output$pca_options <- renderUI({
      ui_pca_eval(d()$cleaned$meta_df, ns = session$ns)
    })
    #-- PCA comparison data for the selected compare mode
    pca_compare_data <- reactive({
      req(input$pca_compare)
      get_pca_compare_data(
        p = p(),
        d = d(),
        pca_compare = input$pca_compare
      )
    })

    #-- Compute PCA once and reuse for both PCA plots

    pca_meta_df <- reactive({
      req(d(), d()$cleaned, d()$cleaned$meta_df)
      d()$cleaned$meta_df
    })
    pca_meta_cols <- reactive({
      req(pca_meta_df())
      unique(c("sample", setdiff(names(pca_meta_df()), "sample")))
    })

    pca_pair_reactive <- reactive({
      req(input$pca_compare)

      pca_p <- p()
      pca_p$pca_compare <- input$pca_compare

      cmp <- pca_compare_data()

      compute_pca_pair(
        before = cmp$before,
        after = cmp$after,
        p = pca_p,
        before_label = "Before",
        after_label = "After",
        meta_cols = pca_meta_cols(),
        meta_df = pca_meta_df(),
        sample_col = "sample"
      )
    })
    #-- PCA score plot
    output$pca_plot <- renderPlot(
      {
        req(input$pca_compare, input$color_col)

        pca_p <- p()
        pca_p$pca_compare <- input$pca_compare
        pca_p$color_col <- input$color_col
        pca_p$shape_col <- input$shape_col

        cmp <- pca_compare_data()

        suppressWarnings({
          plot_pca_from_result(
            p = pca_p,
            pca_pair = pca_pair_reactive(),
            compared_to = cmp$compared_to
          )
        })
      },
      res = 120
    )

    #-- PCA loading plot
    output$pca_loading_plot <- renderPlot(
      {
        req(input$pca_compare)

        cmp <- pca_compare_data()

        plot_pca_loading_from_result(
          pca_pair = pca_pair_reactive(),
          compared_to = cmp$compared_to
        )
      },
      res = 120
    )

    #-- Download all figures as zip folder.
    output$download_fig_zip_btn <- renderUI({
      req(d()$filtered, d()$filtered_corrected)

      download_card(
        "Download Figures",
        "If there are many metabolites, downloading figures may take a few minutes.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "max-width: 250px; display: inline-block;",
            downloadButton(
              outputId = ns("download_fig_zip"),
              label = "Download All Figures",
              class = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })

    # -- progress bar
    progress_reactive <- reactiveVal(0)

    #-- progress for downloading all images
    output$progress_ui <- renderUI({
      req(progress_reactive() > 0, progress_reactive() <= 1)
      div(
        style = "margin-top: 10px;",
        tags$label("Progress:"),
        tags$progress(
          value = progress_reactive(),
          max = 1,
          style = "width: 100%; height: 20px;"
        ),
        tags$span(sprintf("%.0f%%", progress_reactive() * 100))
      )
    })

    output$download_fig_zip <- downloadHandler(
      filename = function() {
        paste0("figures_", Sys.Date(), ".zip")
      },
      content = function(file) {
        .require_pkg("zip", "create a zip archive")

        choices <- list(
          rsd_cal = input$rsd_cal,
          rsd_compare = input$rsd_compare,
          rsd_plot_type = input$rsd_plot_type,
          pca_compare = input$pca_compare,
          color_col = input$color_col,
          fig_format = input$fig_format,
          remove_imputed = p()$remove_imputed,
          transform = p()$transform,
          qcImputeM = p()$qcImputeM,
          samImputeM = p()$samImputeM
        )

        figs <- export_figures(p = choices, d = d(), out_dir = tempdir())

        fig_dir <- normalizePath(figs$fig_dir, winslash = "/", mustWork = TRUE)
        zipfile <- tempfile(fileext = ".zip")
        zip::zipr(zipfile, files = fig_dir)

        file.copy(zipfile, file, overwrite = TRUE)

        unlink(figs$fig_dir, recursive = TRUE, force = TRUE)
        unlink(zipfile, force = TRUE)

        progress_reactive(0)
      }
    )

    #-- Move to next tab after inspecting the corrected data figures
    observeEvent(input$next_export, {
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_export")
    })

    list(
      progress = progress_reactive,
      params = reactive(list(
        rsd_compare = input$rsd_compare,
        rsd_cal = input$rsd_cal,
        rsd_plot_type = input$rsd_plot_type,
        pca_compare = input$pca_compare,
        color_col = input$color_col,
        shape_col = input$shape_col,
        fig_format = input$fig_format
      ))
    )
  })
}
