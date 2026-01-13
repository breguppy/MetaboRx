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
            options = list(container = "body",
                           customClass = "popover-responsive") 
          )
        ),
        uiOutput(ns("met_plot_selectors")),
        width = 400
      ),
      plotOutput(ns("metab_scatter"), height = "600px", width = "600px") %>% withSpinner(color = "#404040"),
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
            shiny::tags$p(shiny::strong("How to read this plot")),
            shiny::tags$p("This page compares relative standard deviation (RSD) in the corrected or transformed and corrected data ",
            "(depending on the setting selected under 'Compare raw data to') to the raw data. ",
            "RSD is computed by dividing the standard deviation of each metabolite by the mean of that metabolite and is expressed ",
            "as a percentage. RSD is computed for each metabolite for QC samples and non-QC samples separtely. RSD can also be computed for non-QC ",
            "samples grouping samples by class type (depending on the settings selected under 'Calculate RSD by')."),
            shiny::tags$hr(),
            shiny::strong("Visualize changes in RSD by: Distrbution"),
            shiny::tags$p(
              "The distributions of RSDs in non-QC samples is displayed in the left panel and the distribution of RSDs in ",
              "QC samples is displayed in the right panel. ",
              "The blue distribution is RSD in the raw data before any correction or transformations is applied. ",
              "The orange distrubution is RSD in the corrected or transformed and corrected data. "
            ),
            shiny::tags$p(
              shiny::strong("Goal: "),
              "after correction/transformation and correction, the orange distributions should be shifted to the left compared to the blue distributions. ",
              "The orange distribution for QC samples should be tall and skinny with the highest density near zero."
            ),
            shiny::tags$hr(),
            shiny::strong("Visualize changes in RSD by: Scatter Plot"),
            shiny::tags$p(
              "In the scatter plot comparison the x-axis is RSD before correction/transformation and correction and the y-axis is RSD after. ",
              "RSDs for non-QC samples are displayed in the left panel and QC samples in the right panel. ",
              "Red dots indicate that RSD increased after correction/transformation and correction. ", 
              "Gray dot indicate no change in RSD after correction/transformation and correction. ",
              "Green dots indicate a decrease in RSD after correction/transformation. ",
              "The percentages of increased, no change, and decreased RSDs are shown at the top of each panel."
            ),
            shiny::tags$p(
              shiny::strong("Goal: "),
              "after correction/transofrmation and correction, the majority of RSDs should decrease for QC samples. ",
              "Non-QC sample RSDs may or may not decrease dramatically after correction/transformation and correction."
            ),
            title = "What information does the RSD comparison plots show?",
            placement = "auto",
            options = list(container = "body",
                           customClass = "popover-responsive") 
          )
        ),
        ui_rsd_eval(ns),
        shiny::tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
          shiny::tags$strong("Metric guide"),
          bslib::popover(
            shiny::tags$button(
              type = "button",
              class = "btn btn-link p-0",
              style = "text-decoration:none;",
              shiny::icon("circle-info")
            ),
            shiny::tags$p(shiny::strong("")),
            shiny::tags$p("The following table show the average and median change in (\u0394) RSD for both QC samples and non-QC samples.",
                          "We include median as a more robust measure of \u0394 RSD."),
            shiny::tags$p("\u0394 RSD = After RSD - Before RSD. "),
            shiny::tags$p(
              shiny::strong("Goal: "),
              "after correction/transformation and correction, RSD should decrease for both QC and non-QC samples. ",
              " In this situation, a more negative number is disirable for all four \u0394 metrics."
            ),
            title = "What metrics are used to evaluate RSD?",
            placement = "auto",
            options = list(container = "body",
                           customClass = "popover-responsive") 
          )
        ),
        uiOutput(ns("rsd_comparison_stats")),
        width = 400
      ),
      plotOutput(ns("rsd_comparison_plot"), height = "540px", width = "900px") %>% withSpinner(color = "#404040")
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
            shiny::tags$p(shiny::strong("What is principal component analysis (PCA)?")),
            shiny::tags$p("PCA is a dimension reduction technique that projects the original data onto components that capture the maxium variance in the data. ",
                          "Principal conponent 1 (PC1) represents the most variance in the data. After PC1, PC2 represents the most variance in the remaining ",
                          "data."),
            shiny::tags$hr(),
            shiny::strong("PCA score plots"),
            shiny::tags$p(
              "The left panel is the 2D PC plot for the raw data and the right panel is the 2D PC plot for the corrected/transformed and corrected data. ",
              "The x-axis is PC1 and y-axis is PC2. The percentage in the parentheses on the axis labels is the variance explained for each conponent. ",
              "Dots in this figure represent samples."
            ),
            shiny::tags$p(
              shiny::strong("Goal: "),
              "after correction/transformation and correction, biological variation should dominate technical variation and signal drift should ",
              "not be visible in right panel. "
              ),
            shiny::tags$ul(
              shiny::tags$li(shiny::strong("When coloring the plot by class: "), "QC samples should cluster together in the right panel."),
              shiny::tags$li(shiny::strong("When coloring the plot by batch or order: "), "there should be no distinct color patterns in the right panel if samples were run using a random injection ordering")
            ),
            shiny::tags$hr(),
            shiny::strong("PCA loading plots"),
            shiny::tags$p(
              "The loading values show how much a metabolite contributes to that PC and the top 10 metabolites for each PC are shown below the PCA plot. ",
              "The magnitude of the loading corresponds to the metabolite's strength of correlation to that PC. ",
              "A metabolite with a large magnitude (close to 1 or -1) has a strong influence/contribution to that PC ",
              "and a metabolite with a small magnitude close to 0 has weak influence/contribution to that PC. ",
              "A positive loading (green) means that a high value in that metabolite corresponds to a high value in that PC. ",
              "A negative loading (red) means a high value in that metabolite corresponds to a low value in that PC."
            ),
            title = "What information does the PCA plots and loading plots show?",
            placement = "auto",
            options = list(container = "body",
                           customClass = "popover-responsive") 
          )
        ),
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
            help = c("If there are many metabolites, downloading figures may take a few minutes."),
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
           rsd_plot_type = input$rsd_plot_type,
           pca_compare = input$pca_compare,
           color_col   = input$color_col,
           fig_format  = input$fig_format
        ))
    )
  })
}