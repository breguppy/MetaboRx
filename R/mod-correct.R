#' Correction module
#'
#' @keywords internal
#' @noRd

mod_correct_ui <- function(id) { 
  ns <- NS(id); 
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
        column(3, tags$h5("Choose Correction Method"), uiOutput(ns("correctionMethod"))),
        column(3, tags$h5("Unavailable Options"), uiOutput(ns("unavailable_options"))),
        actionButton(ns("correct"), "Correct Data with Selected Settings",
                     class="btn-primary btn-lg", width="100%"),
        div(style="margin:12px 0 0 0;", withSpinner(uiOutput(ns("cor_spinner")),
                                                    color="#404040", size=0.6, proxy.height="22px"))
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.2 Post-Correction Filtering",
          uiOutput(ns("post_cor_filter_block")),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Corrected RSD Summary",
            uiOutput(ns("download_cor_rsd_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
            help = c("Creates Excel file with RSD summary before and after correction for samples and QCs. "),
            width = 400,
            position = "right"
          ),
        uiOutput(ns("post_cor_filter_info")) %>% withSpinner(color = "#404040") 
        )
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.3 Post-Correction Transformation",
          uiOutput(ns("transform_block")),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
             title = "Download Transformed RSD Summary",
          uiOutput(ns("download_tc_rsd_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
          help = c("Creates Excel file with RSD summary before and after correction and transformation for samples and QCs."),
          width = 400,
          position = "right"
          ),
         ui_table_scroll("cor_data", ns) %>% withSpinner(color = "#404040")
        )
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.4 Candidate Extreme Values",
          shiny::tags$div(
            style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
            shiny::tags$strong("How detection works"),
            bslib::popover(
              shiny::tags$button(
                type = "button",
                class = "btn btn-link p-0",
                style = "text-decoration:none;",
                shiny::icon("circle-info")
              ),
              shiny::tagList(
                shiny::tags$p(
                  "This screen flags potential extreme values using a 2D PCA / Hotelling T² approach fit on non-QC samples.",
                  "Hotelling’s T² here is computed as the squared Mahalanobis distance in PC1–PC2 space using a PCA model fit on non-QC samples:"
                ),
                shiny::tags$ol(
                  shiny::tags$li(
                    shiny::strong("Log and Scale Metabolites for non-QC samples: "),
                    "Applies log2(x + 1), then standardizes using pooled non-QC samples only."
                  ),
                  shiny::tags$li(
                    shiny::strong("PCA fit (non-QC only): "),
                    "Fits PCA on pooled non-QC rows with complete metabolite data; uses PC1–PC2."
                  ),
                  shiny::tags$li(
                    shiny::strong("T² in PC space for all samples: "),
                    "Projects all complete rows (QC + non-QC) into PC1–PC2 and computes a squared Mahalanobis distance (Hotelling T²)."
                  ),
                  shiny::tags$li(
                    shiny::strong("Ellipse cutoff: "),
                    "Flags samples outside the (1 − α) ellipse using a χ² cutoff with df = 2 (default α = 0.05 → 95%)."
                  ),
                  shiny::tags$li(
                    shiny::strong("Dual z-score rule (only for outlier samples): "),
                    "Within samples outside the ellipse, flags metabolite values only when BOTH ",
                    shiny::strong("|global z| ≥ 3"),
                    " (pooled non-QC scaling) AND ",
                    shiny::strong("|class z| ≥ 3"),
                    " (within that non-QC class)."
                  )
                ),
                shiny::tags$hr(),
                shiny::tags$p(
                  shiny::strong("Interpretation: "),
                  "Red points are samples outside the ellipse. The table reports the specific metabolite values that also satisfy the dual z-score threshold."
                ),
                shiny::tags$p(shiny::strong("Caution: "),
                              "Candidate extreme values are displayed for the user's benefit. ",
                              "Further investigation and justification is needed before categorizing an extreme value as an outlier and removing it.")
              ),
              title = "Candidate extreme value detection",
              placement = "auto",
              options = list(container = "body", customClass = "popover-responsive")
            )
          ),
          
          ui_detect_outliers_options(ns),
          width = 400
      ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Extreme Value Summary",
            uiOutput(ns("download_ev_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
            help = c("Creates Excel file with summary of extreme value detection."),
            width = 400,
            position = "right"
          ),
          uiOutput(ns("outliers_table"))
        )
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.5 Post-Correction/Transformation Metabolite Correlation",
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
              shiny::tags$p("To investigate linear relationships between metabolites, Pearson's r is computed for each pair. A strong positive linear correlation (Pearson's r near 1) means that as one metabolite increases, the other metabolite consistently increases proportionally."),
              shiny::tags$p("All pairwise correlations are computed, but we only allow pairs with a strong positive linear correlations to be displayed here."),
              shiny::tags$p("To view all pairwise correlations, download the Excel displayed on the right."),
              title = "Pearson's r correlations",
              placement = "auto",
              options = list(container = "body",
                             customClass = "popover-responsive") 
            )
          ),
          ui_tc_corr_slider(ns),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Corrected/Transformed Data Metabolite Correlations",
            uiOutput(ns("download_tc_corr_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
            help = c("Creates Excel file with all pairwise metabolite correlations in the raw data and corrected/transformed data."),
            width = 400,
            position = "right"),
          uiOutput(ns("compute_tc_corr_ui")),
          div(style="margin:12px 0 0 0;", withSpinner(uiOutput(ns("tc_corr_spinner")),
                                                      color="#404040")),
          uiOutput(ns("tc_corr_range_info"))
        )
        )
      ),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "2.6 Identify Control Group",
        help = c(
          "If your data contains a control group please select the name of the control group in the dropdown menu on the right. Fold changes will be computed and added to a separate tab in the corrected data Excel file.",
          "If you data does NOT contain a control group check the 'No control group' box."
        )
      ),
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "Download Corrected and Transformed Data",
          tooltip(
            checkboxInput(
              ns("keep_corrected_qcs"),
              "Include QCs in corrected data file",
              FALSE
            ),
            "Check the box if you want corrected QC values in the downloaded corrected data file.",
            placement = "right"
          ),
          uiOutput(
            ns("download_corr_btn"),
            container = div,
            style = "position: absolute; bottom: 15px; right: 15px;"
          ),
          help = c("Creates Excel file with correction settings, corrected data, transformed data, group statistics, fold changes, and MetaboAnalyst Ready tabs."),
          width = 400,
          position = "right"
        ),
            tooltip(
              checkboxInput(ns("no_control"), "No control group", FALSE),
              "Check the box if The data does not have a control group.",
              placement = "right"
            ),
            conditionalPanel(condition = sprintf("!input['%s']", ns("no_control")), uiOutput(ns(
              "control_class_selector"
            )))
         
      )
    )), 
    card(actionButton(ns("next_visualization"), "Next: Evaluate and Visualize Correction",
                      class="btn-primary btn-lg"))
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
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_cor_rsd_data"),
            label    = "Download Corrected RSD Summary",
            class    = "btn btn-secondary"
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
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_tc_rsd_data"),
            label    = "Download Transformed RSD Summary",
            class    = "btn btn-secondary"
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
    
    #--------- 2.4 Candidate Extreme Values server
    # PCA plot and table with candidate extreme values
    output$outliers_table <- renderUI({
      req(filtered_corrected_r(), transformed_r())
      d <- list(filtered_corrected = filtered_corrected_r(), 
                transformed = transformed_r())
      p <- list(out_data = input$out_data, 
                qcImputeM = input$qcImputeM, 
                samImputeM = input$samImputeM)
      ui_outliers(
        p = p,
        d = d,
        pca_output_id = "hotelling_pca",
        ns = ns
      )
    })
    
    output$hotelling_pca <- shiny::renderPlot({
      req(filtered_corrected_r(), transformed_r())
      p <- list(out_data = input$out_data, 
                qcImputeM = input$qcImputeM, 
                samImputeM = input$samImputeM)
      # Use the same df logic as ui_outliers()
      df <- if (p$out_data == "filtered_cor_data") {
        filtered_corrected_r()$df_no_mv
      } else {
        transformed_r()$df_no_mv
      }
      
      res <- detect_hotelling_nonqc_dual_z(df, p)
      if (!is.null(res$pca_plot)) {
        res$pca_plot
      }
    })
    
    output$download_ev_btn <- renderUI({
      req(transformed_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_ev_data"),
            label    = "Download Extreme Value Summary",
            class    = "btn btn-secondary"
          )
        )
      )
    })
    
    output$download_ev_data <- downloadHandler(
      filename = function() {
        sprintf("extreme_values_%s.xlsx", Sys.Date())
      },
      content = function(file) {
        d <- list(filtered_corrected = filtered_corrected_r(), 
                  transformed = transformed_r())
        p <- list(out_data = input$out_data, 
                  qcImputeM = input$qcImputeM, 
                  samImputeM = input$samImputeM)
        
        outlier_wb <- export_outliers_xlsx(p, d)        
        openxlsx::saveWorkbook(outlier_wb, file, overwrite = TRUE)
      }
    )
    
    #---------- 2.5 Post-Correction/Transformation Correlations
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
      # Always read the key so the UI invalidates when tc_corr_data changes
      key <- tc_corr_key_r()
      req(nzchar(key))
      
      # If we've computed for this key, hide
      if (!is.na(computed_tc_key_r()) && identical(computed_tc_key_r(), key)) {
        return(NULL)
      }
      
      tagList(
        tags$div(
          style = "margin-bottom: 8px; color: #555;",
          "Computing correlations may take a while if the data has many metabolites."
        ),
        actionButton(
          ns("compute_tc_corr"),
          "Compute Metabolite Correlations",
          class = "btn btn-primary btn-lg",
          width = "100%"
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
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_tc_corr_data"),
            label    = "Download Metabolite Correlations",
            class    = "btn btn-secondary"
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
    
    #---------- 2.6 Identify Control Group server
    # require filtered and corrected data ********* This needs updating
    output$control_class_selector <- renderUI({
      req(transformed_r())
      df <- cleaned_r()$df
      classes <- unique(df$class[df$class != "QC"])
      dropdown_choices <- c("Select a class..." = "", classes)
      tooltip(
        selectInput(
          ns("control_class"),
          "Control Group",
          choices = dropdown_choices,
          selected = ""
        ),
        "Name of control samples in class column. This class's average will be used to compute fold changes in the corrected data file.",
        placement = "right"
      )
    })
    
    # button for downloading corrected data.
    output$download_corr_btn <- renderUI({
      req(transformed_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_corr_data"),
            label    = "Download Corrected and Transformed Data",
            class    = "btn btn-secondary"
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
          no_control        = isTRUE(input$no_control),
          control_class     = input$control_class %||% ""
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
    
    #---------- Next: Visualize Data
    observeEvent(input$next_visualization, {
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
      out_data           = input$out_data,
      keep_corrected_qcs = isTRUE(input$keep_corrected_qcs),
      tc_corr_threshold = input$tc_corr_threshold,
      no_control         = isTRUE(input$no_control),
      control_class      = input$control_class %||% ""
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