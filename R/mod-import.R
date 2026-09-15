#' @keywords internal

mod_import_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "1. Import Raw Data",
    value = "tab_import",
    # 1.1 Upload raw data
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "1.1 Upload Raw Data",
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
              report_text_data_req(),
              # shiny::tags$img(
              #  src = "example_data_structure.png",
              #  style = "width: 100%; height: auto; display: block;"
              # ),
              title = "Required data structure and information",
              placement = "auto",
              options = list(
                container = "body",
                customClass = "popover-responsive"
              )
            )
          ),
          ui_file_upload(ns),
          width = 400
        ),
        ui_table_scroll("contents", ns)
      )
    ),
    # 1.2 Raw Data Inspection
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "1.2 Raw Data Inspection",
        shiny::tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
          shiny::tags$strong("Data Inspection"),
          bslib::popover(
            shiny::tags$button(
              type = "button",
              class = "btn btn-link p-0",
              style = "text-decoration:none;",
              shiny::icon("circle-info")
            ),
            report_text_data_inspection(),
            title = "What is cleaned and checked in this section",
            placement = "auto",
            options = list(
              container = "body",
              customClass = "popover-responsive"
            )
          )
        ),
        uiOutput(ns("column_selectors")),
        uiOutput(ns("column_warning")),
        uiOutput(ns("withhold_toggle")),
        uiOutput(ns("n_withhold_ui")),
        uiOutput(ns("withhold_selectors_ui")),
        uiOutput(ns("ui_control_class_selector")),
        width = 400
      ),
      uiOutput(ns("basic_info")) |> withSpinner(color = "#404040")
    )),
    # 1.3 Raw data filtering
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "1.3 Missing Value Filtering",
          #shiny::uiOutput(ns("blank_threshold_controls")),
          shiny::tags$div(
            style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
            shiny::tags$strong("Missing Value Filter"),
            bslib::popover(
              shiny::tags$button(
                type = "button",
                class = "btn btn-link p-0",
                style = "text-decoration:none;",
                shiny::icon("circle-info")
              ),
              report_text_mv_filter(),
              title = "Why filter metabolites based on missing values?",
              placement = "auto",
              options = list(
                container = "body",
                customClass = "popover-responsive"
              )
            )
          ),
          shiny::uiOutput(ns("mv_filter_slider")),
          width = 400
        ),
        #shiny::uiOutput(ns("blank_threshold_info")),
        shiny::uiOutput(ns("filter_info")),
        shiny::uiOutput(ns("download_mv_btn"))
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "1.4 Blank Threshold Filtering",
          shiny::uiOutput(ns("blank_threshold_controls"))
        ),
        shiny::uiOutput(ns("blank_threshold_info")),
        width = 400
      )
    ),
    # Next: Choose Correction Settings
    card(
      uiOutput(ns("next_correction_ui"))
    )
  )
}

mod_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    #---------- 1.1 Upload Raw Data server
    data_raw <- reactive({
      req(input$file1)
      read_raw_data(input$file1$datapath)
    })
    output$contents <- renderTable({
      data_raw()
    })

    #---------- 1.2 Raw Data Inspection server
    # requires raw data to display selection choices
    selections_r <- reactive({
      list(
        sample = input$sample_col %||% "",
        batch  = if (isTRUE(input$single_batch)) "batch" else input$batch_col %||% "",
        class  = input$class_col %||% "",
        order  = input$order_col %||% ""
      )
    }) |> debounce(200)

    output$column_selectors <- renderUI({
      df <- req(data_raw())
      ui_nonmet_cols(names(df), ns = session$ns)
    })

    output$column_warning <- renderUI({
      df <- req(data_raw())
      sel <- selections_r()
      ui_column_warning(
        df,
        c(sel$sample, sel$batch, sel$class, sel$order)
      )
    })

    output$withhold_toggle <- renderUI({
      req(data_raw())
      ui_withhold_toggle(ns = session$ns)
    })

    withheld_ids_r <- reactive({
      if (!isTRUE(input$withhold_cols)) {
        return(character(0))
      }
      n <- input$n_withhold %||% 0
      if (n <= 0) {
        return(character(0))
      }
      paste0("withhold_col_", seq_len(n))
    })

    withheld_r <- reactive({
      ids <- withheld_ids_r()
      if (!length(ids)) {
        return(character(0))
      }
      vals <- vapply(ids, function(id) {
        input[[id]] %||% ""
      }, character(1))
      vals <- unique(vals[nzchar(vals)])
      sel <- selections_r()
      setdiff(vals, c(sel$sample, sel$batch, sel$class, sel$order))
    }) |> debounce(200)

    output$n_withhold_ui <- renderUI({
      df <- req(data_raw())
      max_withhold <- max(ncol(df) - 4, 0)

      if (isTRUE(input$withhold_cols)) {
        numericInput(
          ns("n_withhold"),
          "Number of columns to withhold",
          1,
          1,
          max_withhold
        )
      }
    })

    output$withhold_selectors_ui <- renderUI({
      df <- req(data_raw())
      req(input$n_withhold)
      sel <- selections_r()
      cols <- setdiff(
        names(df),
        c(sel$sample, sel$batch, sel$class, sel$order)
      )
      ids <- withheld_ids_r()
      if (!length(ids)) {
        return(NULL)
      }
      prev_all <- isolate(vapply(ids, function(id) {
        input[[id]] %||% ""
      }, character(1)))
      lapply(seq_along(ids), function(i) {
        id <- ids[i]
        prev <- prev_all[i]
        other <- setdiff(prev_all, prev)
        choices_i <- c("Select a column..." = "", setdiff(cols, other))
        selectInput(
          ns(id),
          paste("Select column to withhold #", i),
          choices = choices_i,
          selected = if (nzchar(prev) && prev %in% choices_i) {
            prev
          } else {
            ""
          }
        )
      })
    })

    cleaned_r <- reactive({
      df <- req(data_raw())
      sel <- selections_r()
      col_warn <- ui_column_warning(
        df,
        c(sel$sample, sel$batch, sel$class, sel$order)
      )
      req(is.null(col_warn))
      withheld <- withheld_r()
      # column warings must be NULL for cleaning to happen.
      req(all(nzchar(
        c(sel$sample, sel$batch, sel$class, sel$order)
      )))
      req(length(unique(
        c(sel$sample, sel$batch, sel$class, sel$order)
      )) == 4)
      clean_data(
        df = df,
        sample = sel$sample,
        batch = sel$batch,
        class = sel$class,
        order = sel$order,
        withheld_cols = withheld
      )
    }) |> bindCache(input$file1$datapath, selections_r(), withheld_r())

    output$basic_info <- renderUI({
      cd <- req(cleaned_r())
      ui_basic_info(cd)
    })
    output$ui_control_class_selector <- renderUI({
      cd <- req(cleaned_r())
      ui_control_class_selector(cd$df, ns = session$ns)
    })

    #---------- 1.3 Raw Data filtering server
    
    # Show blank threshold controls only when blanks/PBs exist.
    output$blank_threshold_controls <- shiny::renderUI({
      cd <- req(cleaned_r())
      
      blank_df <- cd$blank_df
      has_blanks <- !is.null(blank_df) && nrow(blank_df) > 0L
      
      if (!has_blanks) {
        return(NULL)
      }
      
      ui_blank_threshold_controls(
        ns = session$ns,
        threshold = input$blank_threshold %||% 3,
        remove_default = isTRUE(input$remove_blank_threshold_cols)
      )
    })
    
    # Missing-value filter controls.
    output$mv_filter_slider <- shiny::renderUI({
      req(cleaned_r())
      ui_filter_slider(ns = session$ns)
    })
    
    # Combined raw-data filtering:
    # 1) remove metabolites with no detected signal
    # 2) apply missing-value filter
    # 3) detect blank-threshold failures among retained metabolites
    # 4) optionally remove metabolites below the blank threshold
    filtered_r <- shiny::reactive({
      cd <- req(cleaned_r())
      
      mv_cutoff <- req(input$mv_cutoff)
      qc_mv_cutoff <- req(input$qc_mv_cutoff)
      filter_rule <- req(input$mv_filter_rule)
      
      # ---------------------------------------------------------------------------
      # 1. Remove metabolites with no detected signal
      # ---------------------------------------------------------------------------
      
      # Metabolites that are all missing/zero in both study and QC samples.
      not_detected_mets <- intersect(
        cd$all_missing_zero_non_qc_cols,
        cd$all_missing_zero_qc_cols
      )
      
      # Metabolites all missing/zero in only one sample type.
      all_missing_zero_non_qc_cols <- setdiff(
        cd$all_missing_zero_non_qc_cols,
        not_detected_mets
      )
      
      all_missing_zero_qc_cols <- setdiff(
        cd$all_missing_zero_qc_cols,
        not_detected_mets
      )
      
      # Remove all metabolites that are completely absent from study samples,
      # QC samples, or both before applying percentage-based filtering.
      drops <- union(
        cd$all_missing_zero_non_qc_cols,
        cd$all_missing_zero_qc_cols
      )
      
      df_for_filtering <- cd$df[
        ,
        !(names(cd$df) %in% drops),
        drop = FALSE
      ]
      
      # ---------------------------------------------------------------------------
      # 2. Missing-value filtering
      # ---------------------------------------------------------------------------
      
      metab_cols <- setdiff(
        names(df_for_filtering),
        c("sample", "batch", "class", "order")
      )
      
      mv_filter_result <- filter_by_missing(
        df = df_for_filtering,
        metab_cols = metab_cols,
        mv_cutoff = mv_cutoff,
        qc_mv_cutoff = qc_mv_cutoff,
        filter_rule = filter_rule
      )
      
      # At this point, these are the metabolites that survived missing-value
      # filtering.
      df_after_missing <- mv_filter_result$df
      
      retained_metab_cols <- setdiff(
        names(df_after_missing),
        c("sample", "batch", "class", "order")
      )
      
      # ---------------------------------------------------------------------------
      # 3. Blank-threshold detection
      # ---------------------------------------------------------------------------
      
      blank_df <- cd$blank_df
      
      has_blanks <-
        !is.null(blank_df) &&
        nrow(blank_df) > 0L
      
      blank_threshold_result <- NULL
      removed_blank_threshold_cols <- character(0)
      
      if (
        has_blanks &&
        length(retained_metab_cols) > 0L
      ) {
        blank_threshold_result <- detect_blank_threshold(
          df = df_after_missing,
          blank_df = blank_df,
          metab_cols = retained_metab_cols,
          class_col = "class",
          qc_label = "QC",
          threshold = input$blank_threshold %||% 3,
          internal_standard_pattern = "ISTD|ITSD"
        )
      }
      
      # ---------------------------------------------------------------------------
      # 4. Optional blank-threshold removal
      # ---------------------------------------------------------------------------
      
      final_df <- df_after_missing
      
      if (
        !is.null(blank_threshold_result) &&
        isTRUE(input$remove_blank_threshold_cols)
      ) {
        blank_filter_result <- apply_blank_threshold_filter(
          df = df_after_missing,
          failed_cols = blank_threshold_result$below_blank_threshold,
          protect_internal_standards = TRUE,
          internal_standard_pattern = "ISTD|ITSD"
        )
        
        final_df <- blank_filter_result$df
        
        removed_blank_threshold_cols <-
          blank_filter_result$removed_blank_threshold_cols
      }
      
      # ---------------------------------------------------------------------------
      # 5. Attach reporting information
      # ---------------------------------------------------------------------------
      
      mv_filter_result$df <- final_df
      
      mv_filter_result$blank_threshold_result <-
        blank_threshold_result
      
      mv_filter_result$removed_blank_threshold_cols <-
        removed_blank_threshold_cols
      
      mv_filter_result$blank_threshold <-
        input$blank_threshold %||% 3
      
      mv_filter_result$remove_blank_threshold_cols <-
        isTRUE(input$remove_blank_threshold_cols)
      
      mv_filter_result$not_detected_mets <-
        not_detected_mets
      
      mv_filter_result$all_missing_zero_non_qc_cols <-
        all_missing_zero_non_qc_cols
      
      mv_filter_result$all_missing_zero_qc_cols <-
        all_missing_zero_qc_cols
      
      mv_filter_result
    })
    
    # Blank-threshold filtering information.
    output$blank_threshold_info <- shiny::renderUI({
      cd <- req(cleaned_r())
      
      blank_df <- cd$blank_df
      has_blanks <- !is.null(blank_df) && nrow(blank_df) > 0L
      
      if (!has_blanks) {
        return(NULL)
      }
      
      fd <- req(filtered_r())
      
      ui_blank_threshold_info(
        blank_threshold_result = fd$blank_threshold_result,
        blank_df = blank_df,
        threshold = fd$blank_threshold,
        remove_blank_threshold_cols = fd$remove_blank_threshold_cols,
        removed_blank_threshold_cols = fd$removed_blank_threshold_cols
      )
    })
    
    # Missing-value filter information.
    output$filter_info <- shiny::renderUI({
      fd <- req(filtered_r())
      ui_filter_info(fd)
    })

    output$download_mv_btn <- renderUI({
      req(filtered_r())
      download_card(
        "Download Missing Value Summary",
        "Creates Excel file with missing value summarized by metabolite, sample, class, batch and class-metabolite.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_mv_data"),
              label    = "Download Missing Value Summary",
              class    = "btn btn-secondary btn-lg"
            )
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
          cleaned = req(cleaned_r()),
          filtered = req(filtered_r())
        )

        wb <- export_mv_xlsx(p, d)
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    #---------- Next: Choose Correction Settings server
    # requires raw correlations
    output$next_correction_ui <- renderUI({
      req(filtered_r())
      actionButton(
        ns("next_correction"),
        "Next: Choose Correction Settings",
        class = "btn-primary btn-lg",
        width = "100%"
      )
    })

    observeEvent(input$next_correction, {
      req(filtered_r())
      validate(
        need(!is.null(cleaned_r()), "Missing cleaned data"),
        need(!is.null(filtered_r()), "Missing filtered data")
      )
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_correct")
    })

    #---------- module output
    # Collect all input parameters from this module.
    params_r <- reactive({
      sel <- selections_r()
      list(
        sample_col         = sel$sample,
        batch_col          = sel$batch,
        class_col          = sel$class,
        order_col          = sel$order,
        withheld_cols      = withheld_r(),
        n_withhold         = input$n_withhold %||% 0,
        no_control         = isTRUE(input$no_control),
        control_class      = input$control_class %||% "",
        mv_cutoff          = input$mv_cutoff,
        qc_mv_cutoff       = input$qc_mv_cutoff,
        filter_rule        = input$filter_rule
      )
    })

    list(
      cleaned = cleaned_r,
      filtered = filtered_r,
      params = params_r
    )
  })
}
