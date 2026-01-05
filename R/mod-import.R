#' @keywords internal

mod_import_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "1. Import Raw Data",
    value = "tab_import",
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "1.1 Upload Raw Data",
          ui_file_upload(ns),
          shiny::tags$div(
            style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
            shiny::tags$strong("Data requirements"),
            bslib::popover(
              shiny::tags$button(
                type = "button",
                class = "btn btn-link p-0",
                style = "text-decoration:none;",
                shiny::icon("circle-info")
              ),
              shiny::tags$p("Your upload must match one of these structures:"),
              shiny::tags$p("Columns can be in any order."),
              shiny::tags$img(
                src = image_src <- knitr::image_uri(system.file("www/example_data_structure.png", package = "QCcorrection")),  
                style = "width: 100%; height: auto; display: block;"
              ),
              title = "Required data structure",
              placement = "auto",
              options = list(container = "body",
                             customClass = "popover-responsive") 
            )
          ),
          help = c(
            "Raw data must be on the first sheet of .xls or .xlsx file.",
            "Data must begin and end with QC samples when sorted by injection order."
          ),
          width = 400
        ),
        ui_table_scroll("contents", ns)
      )
    ),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "1.2 Select Non-metabolite Columns",
        uiOutput(ns("column_selectors")),
        uiOutput(ns("column_warning")),
        ui_withhold_toggle(ns),
        uiOutput(ns("n_withhold_ui")),
        uiOutput(ns("withhold_selectors_ui")),
        width = 400
      ),
      uiOutput(ns("basic_info"))
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "1.3 Filter Raw Data", 
        ui_filter_slider(ns), 
        width = 400
        ),
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "Download Missing Value Summary", 
          uiOutput(ns("download_mv_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
          help = c("Missing value summary by metabolite, sample, class, and batch.",
                   "Missing value summary can also be downloaded on tab 4. Export All"),
          width = 400,
          position = "right"),
        uiOutput(ns("filter_info"))
      )
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "1.4 Raw Data Metabolite Correlations",
        ui_corr_slider(ns),
        help = c("To investigate linear relationships between metabolites, Pearson's r is computed for each pair.",
                 "All pairwise correlations are computed, but we only allow pairs with a strong positive linear correlation to be displayed here.",
                 "To view all pairwise correlation, download the Excel displayed on the right."),
        width = 400
      ),
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "Download Raw Data Metabolite Correlations",
          uiOutput(ns("download_raw_corr_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
          help = c("Creates Excel file with all pairwise metabolite correlations in the raw data."),
          width = 400,
          position = "right"),
        uiOutput(ns("compute_raw_corr_ui")),
        div(style="margin:12px 0 0 0;", withSpinner(uiOutput(ns("corr_spinner")),
                                                    color="#404040")),
        uiOutput(ns("corr_range_info"))
        )
    )),
    card(
      actionButton(
        ns("next_correction"),
        "Next: Choose Correction Settings",
        class = "btn-primary btn-lg"
      )
    )
  )
}

mod_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    data_raw <- reactive({
      req(input$file1)
      read_raw_data(input$file1$datapath)
    })
    output$contents <- renderTable(data_raw())
    
    selections_r <- reactive({
      list(
        sample = input$sample_col %||% "",
        batch  = if (isTRUE(input$single_batch)) "batch" else input$batch_col %||% "",
        class  = input$class_col  %||% "",
        order  = input$order_col  %||% ""
      )
    }) %>% debounce(200)
    
    output$column_selectors <- renderUI({
      req(data_raw())
      ui_nonmet_cols(names(data_raw()), ns = session$ns)
    })
    output$column_warning <- renderUI({
      req(data_raw())
      sel <- selections_r()
      ui_column_warning(data_raw(),
                        c(sel$sample, sel$batch, sel$class, sel$order))
    })
    
    withheld_ids_r <- reactive({
      if (!isTRUE(input$withhold_cols))
        return(character(0))
      n <- input$n_withhold %||% 0
      if (n <= 0)
        return(character(0))
      paste0("withhold_col_", seq_len(n))
    })
    withheld_r <- reactive({
      ids <- withheld_ids_r()
      if (!length(ids))
        return(character(0))
      vals <- vapply(ids, function(id)
        input[[id]] %||% "", character(1))
      vals <- unique(vals[nzchar(vals)])
      sel <- selections_r()
      setdiff(vals, c(sel$sample, sel$batch, sel$class, sel$order))
    }) %>% debounce(200)
    
    observe({
      req(data_raw())
      max_withhold <- max(ncol(data_raw()) - 4, 0)
      output$n_withhold_ui <- renderUI({
        if (isTRUE(input$withhold_cols))
          numericInput(ns("n_withhold"),
                       "Number of columns to withhold",
                       1,
                       1,
                       max_withhold)
      })
    })
    
    output$withhold_selectors_ui <- renderUI({
      req(data_raw(), input$n_withhold)
      sel <- selections_r()
      cols <- setdiff(names(data_raw()),
                      c(sel$sample, sel$batch, sel$class, sel$order))
      ids <- withheld_ids_r()
      if (!length(ids))
        return(NULL)
      prev_all <- isolate(vapply(ids, function(id)
        input[[id]] %||% "", character(1)))
      lapply(seq_along(ids), function(i) {
        id <- ids[i]
        prev <- prev_all[i]
        other <- setdiff(prev_all, prev)
        choices_i <- c("Select a column..." = "", setdiff(cols, other))
        selectInput(
          ns(id),
          paste("Select column to withhold #", i),
          choices = choices_i,
          selected = if (nzchar(prev) && prev %in% choices_i)
            prev
          else
            ""
        )
      })
    })
    
    cleaned_r <- reactive({
      df  <- req(data_raw())
      sel <- selections_r()
      withheld <- withheld_r()
      req(all(nzchar(
        c(sel$sample, sel$batch, sel$class, sel$order)
      )))
      req(length(unique(
        c(sel$sample, sel$batch, sel$class, sel$order)
      )) == 4)
      clean_data(df, sel$sample, sel$batch, sel$class, sel$order, withheld)
    }) %>% bindCache(reactiveVal(NULL)(), selections_r(), withheld_r())
    
    output$basic_info <- renderUI({
      cd <- cleaned_r()
      req(cd)
      ui_basic_info(cd$df, cd$replacement_counts, cd$non_numeric_cols, cd$duplicate_mets, cd$blank_df, cd$below_blank_threshold)
    })
    
    filtered_r <- reactive({
      cd <- req(cleaned_r())
      filter_by_missing(cd$df, setdiff(names(cd$df), c("sample", "batch", "class", "order")), input$mv_cutoff)
    })
    output$filter_info <- renderUI({
      fd <- filtered_r()
      req(fd)
      ui_filter_info(fd$mv_removed_cols,
                     input$mv_cutoff,
                     fd$qc_missing_mets)
    })
    
    # button for downloading missing value report.
    output$download_mv_btn <- renderUI({
      req(cleaned_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_mv_data"),
            label    = "Download Missing Value Info",
            class    = "btn btn-secondary"
          )
        )
      )
    })
    output$download_mv_data <- downloadHandler(
      filename = function() {
        paste0("missing_value_counts_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        p <- list()
        
        d <- list(
          cleaned = cleaned_r()
        )
        
        wb <- export_mv_xlsx(p, d)
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    # Increment whenever filtered df changes
    filtered_version_r <- reactiveVal(0L)
    
    observeEvent(filtered_r()$df, {
      filtered_version_r(filtered_version_r() + 1L)
    }, ignoreInit = TRUE)
    
    # Store the version we last computed correlations for
    computed_version_r <- reactiveVal(NA_integer_)
    output$compute_raw_corr_ui <- renderUI({
      req(filtered_r())
      v <- filtered_version_r()
      
      if (isTRUE(!is.na(computed_version_r())) && identical(computed_version_r(), v)) {
        return(NULL) # hide after computed, until df changes
      }
      
      tagList(
        tags$div(
          style = "margin-bottom: 8px; color: #555;",
          "Computing correlations may take a while if the data has many metabolites."
        ),
        actionButton(
          ns("compute_raw_corr"),
          "Compute Metabolite Correlations",
          class = "btn-primary btn-lg",
          width = "100%"
        )
      )
    })
    
    raw_correlations_r <- eventReactive(input$compute_raw_corr, {
      df <- isolate(filtered_r()$df)
      metab <- setdiff(names(df), c("sample", "batch", "class", "order"))
      compute_pairwise_metabolite_correlations(df, metab)
    })
    observeEvent(input$compute_raw_corr, ignoreInit = TRUE, {
      shinyjs::disable("compute_raw_corr")
      output$corr_spinner <- renderUI({
        on.exit(shinyjs::enable("compute_raw_corr"), add = TRUE)
        raw_correlations_r(); computed_version_r(filtered_version_r()); NULL
      })
    })
    output$corr_spinner <- renderUI(NULL)
    
    observeEvent(raw_correlations_r(), {
      computed_version_r(filtered_version_r())
    }, ignoreInit = TRUE)
    
    output$corr_range_info <- renderUI({
      all_corr <- req(raw_correlations_r())
      ui_corr_range_info(all_corr, input$corr_threshold)
    })
    output$download_raw_corr_btn <- renderUI({
      req(raw_correlations_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_raw_corr_data"),
            label    = "Download Metabolite Correlations",
            class    = "btn btn-secondary"
          )
        )
      )
    })
    output$download_raw_corr_data <- downloadHandler(
      filename = function() {
        paste0("raw_metabolite_correlations_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        wb <- export_corr_xlsx(raw_correlations_r()) 
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    params_r <- reactive({
      sel <- selections_r()
      list(
        sample_col = sel$sample,
        batch_col = sel$batch,
        class_col  = sel$class,
        order_col = sel$order,
        withheld_cols = withheld_r(),
        n_withhold = input$n_withhold %||% 0,
        mv_cutoff = input$mv_cutoff,
        raw_corr_threshold = input$corr_threshold
      )
    })
    
    observeEvent(input$next_correction, {
      validate(
        need(!is.null(cleaned_r()), "Missing cleaned data"),
        need(!is.null(filtered_r()), "Missing filtered data")
      )
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_correct")
    })
    
    # module output
    list(cleaned  = cleaned_r,
         filtered = filtered_r,
         raw_corr = raw_correlations_r,
         params   = params_r)
  })
}
