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
          ui_post_cor_filter(ns),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Corrected RSD Summary",
            uiOutput(ns("download_cor_rsd_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
            help = c("RSD summary before and after correction for samples and QCs. ",
                     "RSD summary can also be downloaded on tab 4. Export All"),
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
          uiOutput(ns("transform_selection_ui")),
          uiOutput(ns("trn_withhold_ui")),
          uiOutput(ns("trn_withhold_selectors_ui")),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
             title = "Download Transformed RSD Summary",
          uiOutput(ns("download_tc_rsd_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
          help = c("RSD summary before and after correction and transformation for samples and QCs."),
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
          ui_detect_outliers_options(ns),
          help = c("PCA and Hotelling's T^2 95% ellipse are computes in the PC1-PC2 space using only the non-QC samples.",
                   "Samples outside the Hotelling's T^2 ellipse are colored red",
                   "Samples listed in the table are outside the Hotelling's T^2 95% limit AND have at least 1 potential extreme metabolite value meaning global AND class |z| is greater than 3."
        )),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Extreme Value Summary",
            uiOutput(ns("download_ev_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
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
          ui_tc_corr_slider(ns),
          help = c("To investigate linear relationships between metabolites Pearson's r is computed for each pair.",
                   "All pairwise correlations are computed, but we only allow pairs with a strong positive linear correlation to be displayed here.",
                   "To view all pairwise correlation, download the Excel displayed on the right."),
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
          help = c("Create Excel file with correction settings, corrected data, transformed data, fold changes, and MetaboAnalyst ready tabs."),
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
    
    ########## Part 2.1: choose correction settings:
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
    
    
    ########## Part 2.2 Post correction filtering:
    filtered_corrected_r <- reactive({
      req(filtered_r(), corrected_r())
      df_filtered <- filtered_r()$df
      df_corrected <- corrected_r()$df
      
      if (isTRUE(input$post_cor_filter)) {
        fil_cor <- filter_by_qc_rsd(df_filtered, 
                                    df_corrected, 
                                    Inf, 
                                    input$remove_imputed, 
                                    c("sample","batch","class","order"))
      } else {
        fil_cor <- filter_by_qc_rsd(df_filtered, 
                                    df_corrected, 
                                    input$rsd_filter, 
                                    input$remove_imputed, 
                                    c("sample","batch","class","order"))
      }
      
      fil_cor 
    })
    
    output$post_cor_filter_info <- renderUI({
      req(filtered_corrected_r())
      ui_postcor_filter_info(filtered_corrected_r(), 
                             input$remove_imputed, 
                             input$rsd_filter, 
                             input$post_cor_filter)
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
    
    ########### Part 2.3 Post-correction Transformation:
    output$transform_selection_ui <- renderUI({
      req(filtered_corrected_r())
      
      df <- if (isTRUE(input$remove_imputed)) filtered_corrected_r()$df_mv else filtered_corrected_r()$df_no_mv
      mc <- setdiff(names(df), c("sample","batch","class","order"))
      
      ui_post_cor_transform(df, mc, ns = session$ns)
    })
    
    
    transformed_r <- reactive({
      req(filtered_corrected_r())
      
      transform_method <- input$transform %||% "none"
      ex_istd          <- input$ex_ISTD %||% TRUE
      withhold_on      <- isTRUE(input$trn_withhold_checkbox)
      
      if (isTRUE(input$remove_imputed)) {
        df_filtered <- filtered_corrected_r()$df_mv
      } else {
        df_filtered <- filtered_corrected_r()$df_no_mv
      }
      
      withheld <- character(0)
      if (withhold_on && !is.null(input$trn_withhold_n)) {
        for (i in seq_len(input$trn_withhold_n)) {
          col <- input[[paste0("trn_withhold_col_", i)]]
          if (!is.null(col) && nzchar(col) && col %in% names(df_filtered)) {
            withheld <- c(withheld, col)
          }
        }
      }
      
      transform_data(filtered_corrected_r(), transform_method, withheld, ex_istd)
    })
    
    observe({
      req(corrected_r(), input$trn_withhold_checkbox)
      
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
      req(corrected_r(), input$trn_withhold_n, input$ex_ISTD)
      cols <- names(corrected_r()$df)
      cols <- setdiff(cols, c("sample", "batch", "class", "order"))
      if (input$ex_ISTD) {
        istd <- grep("ISTD", cols, value = TRUE)
        itsd <- grep("ITSD", cols, value = TRUE)
        cols <- setdiff(cols, c(istd, itsd))
      }
      dropdown_choices <- c("Select a column..." = "", cols)
      
      # Generate list of columns to withhold
      lapply(seq_len(input$trn_withhold_n), function(i) {
        selectInput(
          inputId = ns(paste0("trn_withhold_col_", i)),
          label = paste("Select column to withhold #", i),
          choices = dropdown_choices,
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
    
    ######### Part 2.4 Candidate Extreme Values
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
    
    # ---------- Part 2.5: TC correlations button show/hide  ----------
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
    
    # Store which selection we last computed for
    computed_tc_key_r <- reactiveVal(NA_character_)
    
    # Current selection key (changes when user switches tc_corr_data)
    tc_corr_key_r <- reactive({
      paste(
        input$tc_corr_data %||% "filtered_cor_data",
        isTRUE(input$remove_imputed),
        sep = "::"
      )
    })
    
    # Force button to reappear whenever tc_corr_data (or remove_imputed) changes
    observeEvent(input$tc_corr_data, {
      computed_tc_key_r(NA_character_)
    }, ignoreInit = TRUE)
    
    # Also reset when upstream invalidates (helps when data genuinely updates)
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
    
    ############ Part 2.6 Identify Control Group
    output$control_class_selector <- renderUI({
      req(cleaned_r())
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
    
    observeEvent(input$next_visualization, {
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_visualize")
    })
    
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
      no_control         = isTRUE(input$no_control),
      control_class      = input$control_class %||% ""
    ))
    
    list(
      imputed            = imputed_r,
      corrected          = corrected_r,
      filtered_corrected = filtered_corrected_r,
      transformed        = transformed_r,
      params             = correct_params
    )
  })
}