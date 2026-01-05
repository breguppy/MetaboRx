#' @keywords internal


mod_visualize_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "3. Evaluation Metrics and Visualization",
    value = "tab_visualize",
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "3.1 Visualize Correction with Metabolite Scatter Plots",
        uiOutput(ns("met_plot_selectors")),
        width = 400
      ),
      plotOutput(ns("metab_scatter"), height = "600px", width = "600px") %>% withSpinner(color = "#404040"),
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "3.2 RSD Evaluation",
        ui_rsd_eval(ns),
        uiOutput(ns("rsd_comparison_stats")),
        width = 400
      ),
      plotOutput(ns("rsd_comparison_plot"), height = "540px", width = "900px") %>% withSpinner(color = "#404040")
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "3.3 PCA Evaluation",
        ui_pca_eval(ns),
        width = 400
      ),
      plotOutput(ns("pca_plot"), height = "530px", width = "1000px") %>% withSpinner(color = "#404040"),
      plotOutput(ns("pca_loading_plot"), height = "530px", width = "1000px") %>% withSpinner(color = "#404040")
    )),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "3.4 Select Figure Format",
          help = c("Figures will be downloaded in the format selected on the right."),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Figures",
            uiOutput(ns("download_fig_zip_btn")),
            help = c("All figures can also be downloaded on tab 4. Export All"),
            position = "right"
          ),
          ui_fig_format(ns),
          uiOutput(ns("progress_ui"))
        ))
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
      if (isTRUE(p()$remove_imputed)){
        df_cor <- d()$filtered_corrected$df_mv
      } else {
        df_cor <- d()$filtered_corrected$df_no_mv
      }
      raw_cols <- setdiff(names(df_raw), c("sample","batch","class","order"))
      cor_cols <- setdiff(names(df_cor), c("sample","batch","class","order"))
      cols <- intersect(raw_cols, cor_cols)
      validate(need(length(cols) >= 1, "No overlapping metabolites."))
      selectInput(ns("met_col"), "Metabolite column", choices = cols, selected = cols[1])
    })
    
    #-- Metabolite scatter plot
    output$metab_scatter <- renderPlot({
      req(input$met_col)
      make_met_scatter(d(), p(), input$met_col)
    }, res = 120)
    
    #-- RSD comparison plot
    output$rsd_comparison_plot <- renderPlot(execOnResize = FALSE, res = 120,{
      req(input$rsd_compare, input$rsd_cal)
      
      make_rsd_plot(list(rsd_compare = input$rsd_compare, rsd_cal = input$rsd_cal, rsd_plot_type = input$rsd_plot_type, remove_imputed = p()$remove_imputed), d())
    })
    
    output$rsd_comparison_stats <- renderUI({
      req(input$rsd_compare, input$rsd_cal)
      ui_rsd_stats(list(rsd_compare = input$rsd_compare, rsd_cal = input$rsd_cal, remove_imputed = p()$remove_imputed), d())
    })
    
    #-- PCA plot
    output$pca_plot <- renderPlot({
      req(input$pca_compare, input$color_col)
      pca_p <- p()
      pca_p$pca_compare <- input$pca_compare
      pca_p$color_col <- input$color_col
      make_pca_plot(pca_p, d())
    }, res = 120)
    
    output$pca_loading_plot <- renderPlot({
      req(input$pca_compare, input$color_col)
      pca_p <- p()
      pca_p$pca_compare <- input$pca_compare
      pca_p$color_col <- input$color_col
      make_pca_loading_plot(pca_p, d())
    }, res = 120)
    
    #-- Download all figures as zip folder.
    output$download_fig_zip_btn <- renderUI({
      req(d()$transformed)
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_fig_zip"),
            label    = "Download All Figures",
            class    = "btn btn-secondary"
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
          rsd_cal     = input$rsd_cal,
          rsd_compare = input$rsd_compare,
          rsd_plot_type = input$rsd_plot_type,
          pca_compare = input$pca_compare,
          color_col   = input$color_col,
          fig_format  = input$fig_format,
          remove_imputed = p()$remove_imputed
        )
        rv_data <- list(
          filtered           = d()$filtered,
          imputed            = d()$imputed,
          corrected          = d()$corrected,
          filtered_corrected = d()$filtered_corrected,
          transformed        = d()$transformed
        )
        figs <- export_figures(p = choices, d = rv_data, out_dir = tempdir())
        
        fig_dir <- normalizePath(figs$fig_dir, winslash = "/", mustWork = TRUE)
        zipfile <- tempfile(fileext = ".zip")
        zip::zipr(zipfile, files = fig_dir)
        
        file.copy(zipfile, file, overwrite = TRUE)
        
        unlink(figs$fig_dir, recursive = TRUE, force = TRUE)
        unlink(zipfile, force = TRUE)
        
        # Remove progress bar
        progress_reactive(0)
      }
    )
    
    #-- Move to next tab after inspecting the corrected data figures
    observeEvent(input$next_export, {
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_export")
    })
    
    list(progress = progress_reactive, 
         params   = reactive(list(
           rsd_compare = input$rsd_compare,
           rsd_cal     = input$rsd_cal,
           pca_compare = input$pca_compare,
           color_col   = input$color_col,
           fig_format  = input$fig_format
        ))
    )
  })
}