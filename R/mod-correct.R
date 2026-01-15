#' Correction module
#'
#' @keywords internal
#' @noRd

mod_correct_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel(
    title = "2. Correction Settings",
    value = "tab_correct",
    card(
      style = "background-color:#eee;",
      tags$h4("2.1 Choose Correction Settings"),
      uiOutput(ns("qc_missing_value_warning")),
      fluidRow(
        column(3, tags$h5("Impute Missing QC Values"), uiOutput(ns("qcImpute"))),
        column(3, tags$h5("Impute Missing Sample Values"), uiOutput(ns("sampleImpute"))),
        column(3, tags$h5("Choose Correction Method"), uiOutput(ns(
          "correctionMethod"
        ))),
        column(3, tags$h5("Unavailable Options"), uiOutput(ns(
          "unavailable_options"
        ))),
        actionButton(
          ns("correct"),
          "Correct Data with Selected Settings",
          class = "btn-primary btn-lg",
          width = "100%"
        ),
        div(
          style = "margin:12px 0 0 0;",
          withSpinner(
            uiOutput(ns("cor_spinner")),
            color = "#404040",
            size = 0.6,
            proxy.height = "22px"
          )
        )
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(title = "2.2 Post-Correction Filtering", uiOutput(ns(
          "post_cor_filter_block"
        )), width = 400),
        fluidRow(column(
          4, uiOutput(ns("post_cor_filter_info")) %>% withSpinner(color = "#404040")
        ), column(8, uiOutput(
          ns("outliers_table")
        ))),
        fluidRow(column(4, uiOutput(
          ns("download_cor_rsd_btn")
        )), column(8, uiOutput(
          ns("download_ev_btn")
        )))
      )
    ),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(title = "2.3 Post-Correction Transformation", uiOutput(ns(
        "transform_block"
      )), width = 400),
      fluidRow(column(
        8,
        ui_table_scroll("cor_data", ns) %>% withSpinner(color = "#404040"),
      ),
      column(4, uiOutput(ns(
        "download_tc_rsd_btn"
      )), uiOutput(ns(
        "download_corr_btn"
      ))))
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "1.4 Raw Data Metabolite Correlations",
        shiny::tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
          shiny::tags$strong("Pearson's r correlations"),
          bslib::popover(
            shiny::tags$button(
              type = "button",
              class = "btn btn-link p-0",
              style = "text-decoration:none;",
              shiny::icon("circle-info")
            ),
            report_text_correlations(),
            title = "Pearson's r correlations",
            placement = "auto",
            options = list(container = "body",
                           customClass = "popover-responsive") 
          )
        ),
        uiOutput(ns("raw_corr_slider")),
        width = 400
      ),
      fluidRow(
        column(8, 
               uiOutput(ns("compute_raw_corr_ui")),
               div(style="margin:12px 0 0 0;", withSpinner(uiOutput(ns("corr_spinner")),
                                                           color="#404040")),
               uiOutput(ns("corr_range_info"))
        ),
        column(4, 
               uiOutput(ns("download_raw_corr_btn")))
      )
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "2.4 Post-Correction/Transformation Metabolite Correlation",
        shiny::tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
          shiny::tags$strong("Pearson's r correlations"),
          bslib::popover(
            shiny::tags$button(
              type = "button",
              class = "btn btn-link p-0",
              style = "text-decoration:none;",
              shiny::icon("circle-info")
            ),
            report_text_correlations(),
            title = "Pearson's r correlations",
            placement = "auto",
            options = list(container = "body", customClass = "popover-responsive")
          )
        ),
        uiOutput(ns("tc_corr_slider")),
        width = 400
      ),
      fluidRow(column(
        8,
        uiOutput(ns("compute_tc_corr_ui")),
        div(style = "margin:12px 0 0 0;", withSpinner(uiOutput(
          ns("tc_corr_spinner")
        ), color = "#404040")),
        uiOutput(ns("tc_corr_range_info"))
      ), column(4, uiOutput(
        ns("download_tc_corr_btn")
      )))
      
    )),
    card(uiOutput(ns(
      "next_visualization_ui"
    ))
    )
  )}

mod_correct_server <- function(id, data, params) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    d <- reactive(data()) 
    
    filtered_r <- reactive({
      req(d()$filtered)
      d()$filtered
    })
    cleaned_r <- reactive({
      d()$cleaned
    })
    
    #---------- 2.1: Choose Correction Settings server
    # requires filtered data
    output$qc_missing_value_warning <- renderUI({
      df <- filtered_r()$df
      ui_qc_missing_warning(df)
    })
    
    output$qcImpute <- renderUI({
      df <- filtered_r()$df
      mc <- setdiff(names(df), c('sample','batch','class','order'))
      ui_qc_impute(df, mc, ns = session$ns)
    })
    
    output$sampleImpute <- renderUI({
      df <- filtered_r()$df
      mc <- setdiff(names(df), c('sample','batch','class','order'))
      ui_sample_impute(df, mc, ns = session$ns)
    })
    
    output$correctionMethod <- renderUI({
      ui_correction_method(filtered_r()$df, ns = session$ns)
    })
    
    output$unavailable_options <- renderUI({
      df <- filtered_r()$df
      mc <- setdiff(names(df), c('sample','batch','class','order'))
      ui_unavailable_options(df, mc)
    })
    
    metab_cols_r <- reactive({
      setdiff(names(filtered_r()$df), c("sample","batch","class","order"))
    })
    
    has_qc_na_r <- reactive({
      df <- filtered_r()$df; mc <- metab_cols_r()
      any(is.na(dplyr::filter(df, .data$class == "QC")[, mc, drop = FALSE]))
    })
    
    has_sam_na_r <- reactive({
      df <- filtered_r()$df; mc <- metab_cols_r()
      any(is.na(dplyr::filter(df, .data$class != "QC")[, mc, drop = FALSE]))
    })
    
    imputed_r <- reactive({
      df <- filtered_r()$df; mc <- metab_cols_r()

      qc_method  <- input$qcImputeM %||% "nothing_to_impute"
      sam_method <- input$samImputeM %||% "nothing_to_impute"
      
      if (!has_qc_na_r())  qc_method  <- "nothing_to_impute"
      if (!has_sam_na_r()) sam_method <- "nothing_to_impute"
      
      impute_missing(df, mc, qc_method, sam_method)
    })
    
    corrected_r <- eventReactive(input$correct, {
      imputed <- isolate(imputed_r()); mc <- isolate(metab_cols_r())
      correct_data(imputed$df, mc, isolate(input$corMethod))
    })
    
    observeEvent(input$correct, ignoreInit = TRUE, {
      shinyjs::disable("correct")
      output$cor_spinner <- renderUI({
        on.exit(shinyjs::enable("correct"), add = TRUE)
        corrected_r(); NULL
      })
    })
    output$cor_spinner <- renderUI(NULL)
    
    #--------- 2.2 Post-correction filtering server
    # requires corrected data
    output$post_cor_filter_block <- renderUI({
      req(corrected_r())
      tagList(
        ui_post_cor_filter(ns = session$ns),
        NULL
      )
    })
    
    
    filtered_corrected_r <- reactive({
      req(filtered_r(), corrected_r())
      
      df_filtered  <- filtered_r()$df
      df_corrected <- corrected_r()$df
      
      post_all       <- isTRUE(input$post_cor_filter)          
      remove_imputed <- isTRUE(input$remove_imputed)           
      rsd_cutoff     <- input$rsd_filter %||% Inf              
      
      cutoff_to_use <- if (post_all) Inf else rsd_cutoff
      
      filter_by_qc_rsd(
        df_filtered,
        df_corrected,
        cutoff_to_use,
        remove_imputed,
        c("sample", "batch", "class", "order")
      )
    })
    
    output$post_cor_filter_info <- renderUI({
      req(corrected_r())                # ensures step 2.1 done
      res <- req(filtered_corrected_r())
      
      remove_imputed <- isTRUE(input$remove_imputed)
      rsd_filter     <- input$rsd_filter %||% Inf
      post_cor_all   <- isTRUE(input$post_cor_filter)
      
      ui_postcor_filter_info(res, remove_imputed, rsd_filter, post_cor_all)
    })
    
    
    output$download_cor_rsd_btn <- renderUI({
      req(filtered_corrected_r())
      download_card(
        "Download RSD Summary",
        "Creates Excel file with RSDs of both raw and corrected data for both samples and QCs.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_cor_rsd_data"),
              label    = "Download Corrected RSD Summary",
              class    = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })
    
    output$download_cor_rsd_data <- downloadHandler(
      filename = function() {
        sprintf("corrected_rsd_stats_%s.xlsx", Sys.Date())
      },
      content = function(file) {
        p <- list(
          rsd_compare = "filtered_cor_data",
          remove_imputed = input$remove_imputed
        )
        
        d <- list(
          filtered_corrected = filtered_corrected_r(),
          filtered           = filtered_r()
        )
        
        stats_wb <- export_stats_xlsx(p, d)
        openxlsx::saveWorkbook(stats_wb, file, overwrite = TRUE)
      }
    )
    output$outliers_table <- renderUI({
      req(filtered_corrected_r())
      d <- list(filtered_corrected = filtered_corrected_r())
      p <- list(qcImputeM = input$qcImputeM, 
                samImputeM = input$samImputeM)
      ui_outliers(
        p = p,
        d = d,
        pca_output_id = "hotelling_pca",
        ns = ns
      )
    })
    
    output$hotelling_pca <- shiny::renderPlot({
      req(filtered_corrected_r())
      p <- list(qcImputeM = input$qcImputeM, 
                samImputeM = input$samImputeM)
      df <- filtered_corrected_r()$df_no_mv
      
      res <- detect_hotelling_nonqc_dual_z(df, p)
      if (!is.null(res$pca_plot)) {
        res$pca_plot
      }
    })
    
    output$download_ev_btn <- renderUI({
      req(filtered_corrected_r())
      download_card(
        "Download Extreme Value Summary",
        "Creates Excel file with summary of extreme value detection.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_ev_data"),
              label    = "Download Extreme Value Summary",
              class    = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })
    
    output$download_ev_data <- downloadHandler(
      filename = function() {
        sprintf("extreme_values_%s.xlsx", Sys.Date())
      },
      content = function(file) {
        d <- list(filtered_corrected = filtered_corrected_r())
        p <- list(qcImputeM = input$qcImputeM, 
                  samImputeM = input$samImputeM)
        
        outlier_wb <- export_outliers_xlsx(p, d)        
        openxlsx::saveWorkbook(outlier_wb, file, overwrite = TRUE)
      }
    )
    
    #---------- 2.3 Post-correction Transformation server
    # Requires filtered and corrected data
    output$transform_block <- renderUI({
      req(filtered_corrected_r())
      tagList(
        uiOutput(ns("transform_selection_ui")),
        uiOutput(ns("trn_withhold_ui")),
        uiOutput(ns("trn_withhold_selectors_ui"))
      )
    })
    
    output$transform_selection_ui <- renderUI({
      req(filtered_corrected_r())
      
      df <- if (isTRUE(input$remove_imputed)) filtered_corrected_r()$df_mv else filtered_corrected_r()$df_no_mv
      mc <- setdiff(names(df), c("sample","batch","class","order"))
      
      ui_post_cor_transform(df, mc, ns = session$ns)
    })
    
    transformed_r <- reactive({
      req(filtered_corrected_r())
      
      transform_method <- input$transform %||% "none"
      ex_istd          <- isTRUE(input$ex_ISTD)
      withhold_on      <- isTRUE(input$trn_withhold_checkbox)
      
      df_filtered <- if (isTRUE(input$remove_imputed)) {
        filtered_corrected_r()$df_mv
      } else {
        filtered_corrected_r()$df_no_mv
      }
      
      withheld <- character(0)
      n_withhold <- input$trn_withhold_n %||% 0L
      
      if (withhold_on && n_withhold > 0L) {
        for (i in seq_len(n_withhold)) {
          col <- input[[paste0("trn_withhold_col_", i)]] %||% ""
          if (nzchar(col) && col %in% names(df_filtered)) {
            withheld <- c(withheld, col)
          }
        }
      }
      
      transform_data(filtered_corrected_r(), transform_method, withheld, ex_istd)
    })
    
    
    observe({
      req(filtered_corrected_r(), input$trn_withhold_checkbox)
      
      max_withhold <- max(ncol(corrected_r()$df) - 4, 0)
      
      output$trn_withhold_ui <- renderUI({
        if (input$transform == "TRN") {
          numericInput(
            inputId = ns("trn_withhold_n"),
            label = "Number of columns to withold from TRN",
            value = 1,
            min = 1,
            max = max_withhold
          )
        }
      })
    })
    
    output$trn_withhold_selectors_ui <- renderUI({
      req(filtered_corrected_r())
      n_withhold <- input$trn_withhold_n %||% 0L
      if (n_withhold <= 0L || !identical(input$transform %||% "none", "TRN")) return(NULL)
      
      ex_istd <- isTRUE(input$ex_ISTD)
      
      cols <- setdiff(names(corrected_r()$df), c("sample","batch","class","order"))
      if (ex_istd) {
        cols <- setdiff(cols, c(grep("ISTD", cols, value = TRUE), grep("ITSD", cols, value = TRUE)))
      }
      
      dropdown_choices <- c("Select a column..." = "", cols)
      
      lapply(seq_len(n_withhold), function(i) {
        selectInput(
          inputId  = ns(paste0("trn_withhold_col_", i)),
          label    = paste("Select column to withhold #", i),
          choices  = dropdown_choices,
          selected = ""
        )
      })
    })
    
    
    output$cor_data <- renderTable({
      req(transformed_r())
      if (isTRUE(input$remove_imputed)) {
        df <- transformed_r()$df_mv
      } else {
        df <- transformed_r()$df_no_mv
      }
      df
    })
    
    output$download_tc_rsd_btn <- renderUI({
      req(transformed_r())
      download_card(
        "Download Transformed RSD Summary",
        "Creates Excel file with RSD summary before and after correction and transformation for samples and QCs.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_tc_rsd_data"),
              label    = "Download Transformed RSD Summary",
              class    = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })
    
    output$download_tc_rsd_data <- downloadHandler(
      filename = function() {
        sprintf("transformed_rsd_stats_%s.xlsx", Sys.Date())
      },
      content = function(file) {
        p <- list(
          rsd_compare = "transformed_cor_data",
          remove_imputed = input$remove_imputed
        )
        
        d <- list(
          filtered_corrected = filtered_corrected_r(),
          filtered           = filtered_r(),
          transformed        = transformed_r()
        )
        
        stats_wb <- export_stats_xlsx(p, d)                
        openxlsx::saveWorkbook(stats_wb, file, overwrite = TRUE)
      }
    )
    
    output$download_corr_btn <- renderUI({
      req(transformed_r())
      download_card(
        "Download Corrected and Transformed Data",
        htmltools::tagList(
          tooltip(
          checkboxInput(
            ns("keep_corrected_qcs"),
            "Include QCs in corrected data file",
            FALSE
          ),
          "Check the box if you want corrected QC values in the downloaded corrected data file.",
          placement = "right"
        ),
        htmltools::tags$p("Creates Excel file with correction settings, corrected data, ",
                          "transformed data, group statistics, fold changes, and MetaboAnalyst Ready tabs.")
        ),
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_corr_data"),
              label    = "Download Corrected and Transformed Data",
              class    = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })
    
    output$download_corr_data <- downloadHandler(
      filename = function() {
        paste0("corrected_data_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        fc <- isolate(filtered_corrected_r())
        tr <- isolate(transformed_r())
        cr <- isolate(corrected_r())
        p_in <- params()  
        
        p <- list(
          sample_col        = p_in$sample_col,
          batch_col         = p_in$batch_col,
          class_col         = p_in$class_col,
          order_col         = p_in$order_col,
          Frule             = p_in$Frule,
          remove_imputed    = isTRUE(input$remove_imputed),
          rsd_cutoff        = fc$rsd_cutoff,
          transform         = input$transform,
          ex_ISTD           = isTRUE(input$ex_ISTD),
          keep_corrected_qcs= isTRUE(input$keep_corrected_qcs),
          tc_corr_threshold = input$tc_corr_threshold,
          no_control        = isTRUE(p_in$no_control),
          control_class     = p_in$control_class
        )
        
        rv <- list(
          cleaned            = cleaned_r(),
          filtered           = filtered_r(),
          imputed            = imputed_r(),
          corrected          = cr,
          filtered_corrected = fc,
          transformed        = tr
        )
        
        wb <- export_xlsx(p, rv)                      
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    #---------- 2.4 Post-Correction/Transformation Correlations
    #---------- 1.4 Raw Data Metabolite Correlations server
    # requires missing value filtered data
    
    # move this to mod_correct so all correlations can be exported together
    output$raw_corr_slider <- renderUI({
      req(filtered_r())
      ui_raw_corr_slider(ns = session$ns)
    })
    
    output$compute_raw_corr_ui <- renderUI({
      req(filtered_r())
      v <- filtered_version_r()
      
      if (isTRUE(!is.na(computed_version_r())) && identical(computed_version_r(), v)) {
        return(NULL) # hide after computed, until df changes
      }
      
      tagList(
        tags$div(
          style = "width: 100%; text-align: center;",
          tags$div(
            style = "max-width: 350px; display: inline-block;",
            actionButton(
              ns("compute_raw_corr"),
              "Compute Metabolite Correlations",
              class = "btn-primary btn-lg",
              width = "100%"
            )
          )
        ),
        tags$div(
          style = "margin-bottom: 8px; color: #555;",
          "Computing correlations may take a while if the data has many metabolites."
        ),
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
      download_card(
        "Download Raw Data Metabolite Correlations",
        "Creates Excel file with all pairwise metabolite correlations in the raw data.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_raw_corr_data"),
              label    = "Download Metabolite Correlations",
              class    = "btn btn-secondary btn-lg"
            )
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
    output$tc_corr_slider <- renderUI({
      req(transformed_r())
      ui_tc_corr_slider(ns = session$ns)
    })
    tc_corr_input_df_r <- reactive({
      req(filtered_corrected_r(), transformed_r())
      
      use_mv <- isTRUE(input$remove_imputed)
      use_filtered <- identical(input$tc_corr_data, "filtered_cor_data")
      
      if (use_filtered) {
        if (use_mv) filtered_corrected_r()$df_mv else filtered_corrected_r()$df_no_mv
      } else {
        if (use_mv) transformed_r()$df_mv else transformed_r()$df_no_mv
      }
    })
    
    computed_tc_key_r <- reactiveVal(NA_character_)
    
    tc_corr_key_r <- reactive({
      paste(
        input$tc_corr_data %||% "filtered_cor_data",
        isTRUE(input$remove_imputed),
        sep = "::"
      )
    })
    
    observeEvent(input$tc_corr_data, {
      computed_tc_key_r(NA_character_)
    }, ignoreInit = TRUE)
    
    observeEvent(list(filtered_corrected_r(), transformed_r()), {
      computed_tc_key_r(NA_character_)
    }, ignoreInit = TRUE)
    
    output$compute_tc_corr_ui <- renderUI({
      req(transformed_r())
      # Always read the key so the UI invalidates when tc_corr_data changes
      key <- tc_corr_key_r()
      req(nzchar(key))
      
      # If we've computed for this key, hide
      if (!is.na(computed_tc_key_r()) &&
          identical(computed_tc_key_r(), key)) {
        return(NULL)
      }
      
      tagList(
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "max-width: 350px; display: inline-block;",
            actionButton(
              ns("compute_tc_corr"),
              "Compute Metabolite Correlations",
              class = "btn btn-primary btn-lg",
              width = "100%"
            )
          )
        ),
        
        tags$div(
          style = "margin-bottom: 8px; color: #555;",
          "Computing correlations may take a while if the data has many metabolites."
        )
      )
    })
    
    tc_correlations_r <- eventReactive(input$compute_tc_corr, {
      df <- req(tc_corr_input_df_r())
      metab <- setdiff(names(df), c("sample", "batch", "class", "order"))
      compute_pairwise_metabolite_correlations(df, metab)
    })
    
    observeEvent(input$compute_tc_corr, ignoreInit = TRUE, {
      shinyjs::disable(ns("compute_tc_corr"))
      output$tc_corr_spinner <- renderUI({
        on.exit(shinyjs::enable(ns("compute_tc_corr")), add = TRUE)
        
        tc_correlations_r()
        
        computed_tc_key_r(tc_corr_key_r())
        NULL
      })
    })
    
    output$tc_corr_spinner <- renderUI(NULL)
    
    output$tc_corr_range_info <- renderUI({
      all_corr <- req(tc_correlations_r())
      ui_corr_range_info(all_corr, input$tc_corr_threshold)                
    })
    
    output$download_tc_corr_btn <- renderUI({
      req(tc_correlations_r())
      download_card(
        "Download Corrected/Transformed Data Metabolite Correlations",
        "Creates Excel file with all pairwise metabolite correlations in the raw data and corrected/transformed data.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_tc_corr_data"),
              label    = "Download Metabolite Correlations",
              class    = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })
    
    output$download_tc_corr_data <- downloadHandler(
      filename = function() {
        paste0("corrected_metabolite_correlations_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        if (input$tc_corr_data == "filtered_cor_data") {
          d_type <- "Corrected"}
        else {
          d_type <- "Transformed and Corrected"
        }
        wb <- export_corr_xlsx(d()$raw_corr, tc_correlations_r(), d_type2 = d_type) 
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    #---------- Next: Visualize Data
    output$next_visualization_ui <- renderUI({
      req(tc_correlations_r()) 
      actionButton(
        ns("next_visualization"), 
        "Next: Evaluate and Visualize Correction",
        class="btn-primary btn-lg"
        )
    })
    observeEvent(input$next_visualization, {
      req(tc_correlations_r())
      validate(
        need(!is.null(filtered_corrected_r()), "Missing corrected data"),
        need(!is.null(transformed_r()), "Missing transformed data data")
      )
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_visualize")
    })
    
    #--------- Module outputs
    correct_params <- reactive(list(
      qcImputeM          = input$qcImputeM %||% "nothing_to_impute",
      samImputeM         = input$samImputeM %||% "nothing_to_impute",
      remove_imputed     = isTRUE(input$remove_imputed),
      post_cor_filter    = input$post_cor_filter,
      rsd_cutoff         = filtered_corrected_r()$rsd_cutoff,
      transform          = input$transform,
      ex_ISTD            = isTRUE(input$ex_ISTD),
      keep_corrected_qcs = isTRUE(input$keep_corrected_qcs),
      tc_corr_threshold = input$tc_corr_threshold
    ))
    
    list(
      imputed            = imputed_r,
      corrected          = corrected_r,
      filtered_corrected = filtered_corrected_r,
      transformed        = transformed_r,
      tc_corr            = tc_correlations_r,
      params             = correct_params
    )
  })
}