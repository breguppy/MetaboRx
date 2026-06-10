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
        column(3, tags$h5("How to choose a correction method"), uiOutput(ns(
          "how_to_correct"
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
        fluidRow(
          column(
            6, uiOutput(ns("post_cor_filter_info")) |> withSpinner(color = "#404040")
          ),
          column(
            6,
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
                report_text_rsd_table(),
                title = "What metrics are used to evaluate RSD?",
                placement = "auto",
                options = list(
                  container = "body",
                  customClass = "popover-responsive"
                )
              )
            ),
            uiOutput(ns("rsd_comparison_stats")),
            uiOutput(ns("download_cor_rsd_btn")),
          )
        )
      )
    ),
    card(
      fluidRow(
        column(
          width = 4,
          htmltools::tags$h4("2.3 Candidate Extreme Values"),
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
              report_text_ev_detection(),
              title = "Candidate extreme value detection",
              placement = "auto",
              options = list(container = "body", customClass = "popover-responsive")
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 7,
          uiOutput(ns("outliers_table"))
        ),
        column(
          width = 5,
          shiny::plotOutput(ns("hotelling_pca"), height = "400px"),
          uiOutput(ns("download_ev_btn"))
        )
      )
    ),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(title = "2.4 Post-Correction Transformation", uiOutput(ns(
        "transform_block"
      )), width = 400),
      fluidRow(
        column(
          8,
          ui_table_scroll("cor_data", ns) |> withSpinner(color = "#404040"),
          uiOutput(ns("download_cor_btn"))
        ),
        column(
          4,
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
              report_text_rsd_tc_table(),
              title = "What metrics are used to evaluate RSD?",
              placement = "auto",
              options = list(
                container = "body",
                customClass = "popover-responsive"
              )
            )
          ),
          uiOutput(ns("post_transform_rsd_compare")),
          uiOutput(ns("download_tc_rsd_btn"))
        )
      ),
    )),
    card(layout_sidebar(
      sidebar = ui_sidebar_block(
        title = "2.5 Metabolite Correlation (Optional)",
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
        uiOutput(ns("correlation_slider")),
        width = 400
      ),
      fluidRow(column(
        8,
        uiOutput(ns("compute_corr_ui")),
        div(style = "margin:12px 0 0 0;", withSpinner(uiOutput(
          ns("corr_spinner")
        ), color = "#404040")),
        uiOutput(ns("corr_range_info"))
      ), column(4, uiOutput(
        ns("download_corr_btn")
      )))
    )),
    card(uiOutput(ns(
      "next_visualization_ui"
    )))
  )
}

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
      mc <- setdiff(names(df), c("sample", "batch", "class", "order"))
      ui_qc_impute(df, mc, ns = session$ns)
    })

    output$sampleImpute <- renderUI({
      df <- filtered_r()$df
      mc <- setdiff(names(df), c("sample", "batch", "class", "order"))
      ui_sample_impute(df, mc, ns = session$ns)
    })

    output$correctionMethod <- renderUI({
      ui_correction_method(filtered_r()$df, ns = session$ns)
    })

    # output$unavailable_options <- renderUI({
    #  df <- filtered_r()$df
    #  ui_unavailable_options(df)
    # })
    output$how_to_correct <- renderUI({
      df <- filtered_r()$df
      ui_how_to_correct(df)
    })

    metab_cols_r <- reactive({
      setdiff(names(filtered_r()$df), c("sample", "batch", "class", "order"))
    })

    has_qc_na_r <- reactive({
      df <- filtered_r()$df
      mc <- metab_cols_r()
      any(is.na(dplyr::filter(df, .data$class == "QC")[, mc, drop = FALSE]))
    })

    has_sam_na_r <- reactive({
      df <- filtered_r()$df
      mc <- metab_cols_r()
      any(is.na(dplyr::filter(df, .data$class != "QC")[, mc, drop = FALSE]))
    })

    imputed_r <- reactive({
      df <- filtered_r()$df
      mc <- metab_cols_r()

      qc_method <- input$qcImputeM %||% "nothing_to_impute"
      sam_method <- input$samImputeM %||% "nothing_to_impute"

      if (!has_qc_na_r()) qc_method <- "nothing_to_impute"
      if (!has_sam_na_r()) sam_method <- "nothing_to_impute"

      impute_missing(df, mc, qc_method, sam_method)
    })

    corrected_r <- eventReactive(input$correct, {
      imputed <- isolate(imputed_r())
      mc <- isolate(metab_cols_r())
      correct_data(imputed$df, mc, isolate(input$corMethod))
    })

    observeEvent(input$correct, ignoreInit = TRUE, {
      shinyjs::disable("correct")
      output$cor_spinner <- renderUI({
        on.exit(shinyjs::enable("correct"), add = TRUE)
        corrected_r()
        NULL
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

      df_filtered <- filtered_r()$df
      df_corrected <- corrected_r()$df

      # remove imputed values if selected:
      if (isTRUE(input$remove_imputed)) {
        df_corrected <- remove_imputed_from_corrected(
          raw_df = df_filtered,
          corrected_df = df_corrected
        )
      }

      pct_threshold <- input$qc_average_pct_threshold %||% 100
      pct_threshold <- as.numeric(pct_threshold)

      if (
        length(pct_threshold) != 1L ||
          is.na(pct_threshold) ||
          pct_threshold < 0
      ) {
        pct_threshold <- 100
      }

      removed_metabolites <- character(0)

      flagged_mets <- get_metabs_pct_diff_vs_qc_average(
        df = df_corrected,
        percent_threshold = pct_threshold
      )

      if (isTRUE(input$remove_qc_average_pct_filter)) {
        average_diff_results <- remove_metabs_pct_diff_vs_qc_average(
          df = df_corrected,
          percent_threshold = pct_threshold,
          return_result = TRUE
        )

        df_corrected <- average_diff_results$df
        removed_metabolites <- average_diff_results$removed_metabolites
      }

      post_all <- isTRUE(input$post_cor_filter)
      remove_imputed <- isTRUE(input$remove_imputed)
      rsd_cutoff <- input$rsd_filter %||% Inf

      cutoff_to_use <- if (post_all) Inf else rsd_cutoff

      filtered_cor_results <- filter_by_qc_rsd(
        raw_df = df_filtered,
        corrected_df = df_corrected,
        rsd_cutoff = cutoff_to_use,
        remove_imputed = remove_imputed,
        metadata_cols = c("sample", "batch", "class", "order")
      )

      filtered_cor_results$percent_threshold <- pct_threshold
      filtered_cor_results$flagged_mets <- flagged_mets
      filtered_cor_results$removed_mets_pct_diff <- removed_metabolites

      filtered_cor_results
    })

    output$post_cor_filter_info <- renderUI({
      req(corrected_r())
      res <- req(filtered_corrected_r())

      remove_imputed <- isTRUE(input$remove_imputed)
      rsd_filter <- input$rsd_filter %||% Inf
      post_cor_all <- isTRUE(input$post_cor_filter)

      ui_postcor_filter_info(res, remove_imputed, rsd_filter, post_cor_all, input$remove_qc_average_pct_filter)
    })

    output$rsd_comparison_stats <- renderUI({
      d <- list(
        filtered_corrected = filtered_corrected_r(),
        filtered = filtered_r()
      )
      ui_rsd_stats(
        compare_to = "filtered_cor_data",
        list(remove_imputed = input$remove_imputed),
        d
      )
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
    hotelling_res_r <- reactive({
      req(filtered_corrected_r())

      df <- filtered_corrected_r()$df_no_mv

      p <- list(
        qcImputeM = input$qcImputeM,
        samImputeM = input$samImputeM
      )

      detect_hotelling_nonqc_dual_z(df, p)
    })

    output$hotelling_pca <- renderPlot({
      res <- req(hotelling_res_r())
      res$pca_plot
    })

    output$outliers_table <- renderUI({
      res <- req(hotelling_res_r())

      ui_outliers_table(
        detect_result = res,
        top_n = 10L,
        sample_col = "sample",
        class_col = "class"
      )
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
        p <- list(
          qcImputeM = input$qcImputeM,
          samImputeM = input$samImputeM
        )

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
      mc <- setdiff(names(df), c("sample", "batch", "class", "order"))

      ui_post_cor_transform(df, mc, ns = session$ns)
    })

    transformed_r <- reactive({
      req(filtered_corrected_r())

      transform_method <- input$transform %||% "none"
      ex_istd <- isTRUE(input$ex_ISTD)
      withhold_on <- isTRUE(input$trn_withhold_checkbox)

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
      if (n_withhold <= 0L || !identical(input$transform %||% "none", "TRN")) {
        return(NULL)
      }

      ex_istd <- isTRUE(input$ex_ISTD)

      cols <- setdiff(names(corrected_r()$df), c("sample", "batch", "class", "order"))
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

    output$post_transform_rsd_compare <- renderUI({
      req(transformed_r())

      d <- list(
        filtered_corrected = filtered_corrected_r(),
        filtered           = filtered_r(),
        transformed        = transformed_r()
      )

      ui_rsd_stats(
        compare_to = "transformed_cor_data",
        list(remove_imputed = input$remove_imputed),
        d
      )
    })

    output$download_tc_rsd_btn <- renderUI({
      req(transformed_r())
      if (!identical(input$transform, "none")) {
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
      } else {
        NULL
      }
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

    output$download_cor_btn <- renderUI({
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
          htmltools::tags$p(
            "Creates Excel file with correction settings, corrected data, ",
            "transformed data, group statistics, fold changes, and MetaboAnalyst Ready tabs."
          )
        ),
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_cor_data"),
              label    = "Download Corrected and Transformed Data",
              class    = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })

    output$download_cor_data <- downloadHandler(
      filename = function() {
        paste0("corrected_data_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        fc <- isolate(filtered_corrected_r())
        tr <- isolate(transformed_r())
        cr <- isolate(corrected_r())
        p_in <- params()

        p <- list(
          sample_col = p_in$sample_col,
          batch_col = p_in$batch_col,
          class_col = p_in$class_col,
          order_col = p_in$order_col,
          Frule = p_in$Frule,
          remove_imputed = isTRUE(input$remove_imputed),
          rsd_cutoff = fc$rsd_cutoff,
          transform = input$transform,
          ex_ISTD = isTRUE(input$ex_ISTD),
          keep_corrected_qcs = isTRUE(input$keep_corrected_qcs),
          tc_corr_threshold = input$tc_corr_threshold,
          no_control = isTRUE(p_in$no_control),
          control_class = p_in$control_class
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

    #---------- 2.4 Metabolite Correlations
    # Requires filtered(), filtered_corrected_r(), and transformed_r()
    output$correlation_slider <- renderUI({
      req(transformed_r())
      ui_correlation_slider(ns = session$ns)
    })

    .pick_df_version <- function(x, remove_imputed) {
      if (isTRUE(remove_imputed)) x$df_mv else x$df_no_mv
    }

    .compute_corr <- function(df) {
      metab <- setdiff(names(df), c("sample", "batch", "class", "order"))
      compute_pairwise_metabolite_correlations(df, metab)
    }

    all_corr_r <- reactiveVal(NULL)

    computed_key_r <- reactiveVal(NA_character_)

    .safe_dim <- function(df) {
      if (is.null(df)) {
        return("NAxNA")
      }
      paste0(nrow(df), "x", ncol(df))
    }

    corr_key_r <- reactive({
      req(filtered_r(), filtered_corrected_r())

      remove_imputed <- isTRUE(input$remove_imputed)
      tr_method <- input$transform %||% "none"
      include_trn <- !identical(tr_method, "none")

      raw_df <- filtered_r()$df
      cor_df <- .pick_df_version(filtered_corrected_r(), remove_imputed)

      # Only touch transformed_r() when it should exist for the key
      trn_dim <- if (include_trn) {
        trn_df <- .pick_df_version(req(transformed_r()), remove_imputed)
        .safe_dim(trn_df)
      } else {
        "SKIP"
      }

      paste0(
        "remove_imputed=", remove_imputed,
        "|transform=", tr_method,
        "|raw_dim=", .safe_dim(raw_df),
        "|cor_dim=", .safe_dim(cor_df),
        "|trn_dim=", trn_dim
      )
    })

    observeEvent(
      list(
        filtered_r(), filtered_corrected_r(), transformed_r(),
        input$remove_imputed, input$transform
      ),
      {
        computed_key_r(NA_character_)
      },
      ignoreInit = TRUE
    )

    output$compute_corr_ui <- renderUI({
      req(filtered_corrected_r())
      key <- corr_key_r()
      req(nzchar(key))

      # IMPORTANT: force dependency on computed_key_r()
      ck <- computed_key_r()

      if (!is.na(ck) && identical(ck, key)) {
        return(NULL)
      }

      tagList(
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "max-width: 350px; display: inline-block;",
            actionButton(
              ns("compute_corr"),
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


    compute_all_correlations_r <- eventReactive(input$compute_corr, {
      req(filtered_r(), filtered_corrected_r(), transformed_r())

      remove_imputed <- isTRUE(input$remove_imputed)
      transform_meth <- input$transform %||% "none"

      raw_df <- filtered_r()$df
      cor_df <- .pick_df_version(filtered_corrected_r(), remove_imputed)

      # Only compute transformed correlations if transform != "none"
      do_transformed <- !identical(transform_meth, "none")
      transformed_df <- if (do_transformed) .pick_df_version(transformed_r(), remove_imputed) else NULL

      list(
        raw = .compute_corr(raw_df),
        corrected = .compute_corr(cor_df),
        transformed = if (do_transformed) .compute_corr(transformed_df) else NULL,
        transformed_included = do_transformed,
        transform_method = transform_meth
      )
    })


    observeEvent(input$compute_corr, ignoreInit = TRUE, {
      shinyjs::disable("compute_corr")
      output$corr_spinner <- renderUI({
        on.exit(shinyjs::enable("compute_corr"), add = TRUE)

        res <- compute_all_correlations_r()
        all_corr_r(res)

        computed_key_r(corr_key_r())
        NULL
      })
    })
    output$corr_spinner <- renderUI(NULL)


    output$corr_range_info <- renderUI({
      all_corr <- req(all_corr_r())
      ui_corr_range_info(all_corr, input$corr_threshold)
    })

    output$download_corr_btn <- renderUI({
      req(all_corr_r())
      download_card(
        "Download Corrected/Transformed Data Metabolite Correlations",
        "Creates Excel file with all pairwise metabolite correlations in the raw data and corrected/transformed data.",
        div(
          style = "width: 100%; text-align: center;",
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = ns("download_corr_data"),
              label    = "Download Metabolite Correlations",
              class    = "btn btn-secondary btn-lg"
            )
          )
        )
      )
    })

    output$download_corr_data <- downloadHandler(
      filename = function() {
        paste0("metabolite_correlations_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        wb <- export_corr_xlsx(compute_all_correlations_r())
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    #---------- Next: Visualize Data
    output$next_visualization_ui <- renderUI({
      # req(all_corr_r())
      actionButton(
        ns("next_visualization"),
        "Next: Evaluate and Visualize Correction",
        class = "btn-primary btn-lg"
      )
    })
    observeEvent(input$next_visualization, {
      # req(all_corr_r())
      validate(
        need(!is.null(filtered_corrected_r()), "Missing corrected data"),
        need(!is.null(transformed_r()), "Missing transformed data data")
      )
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_visualize")
    })

    #--------- Module outputs
    correct_params <- reactive(list(
      qcImputeM = input$qcImputeM %||% "nothing_to_impute",
      samImputeM = input$samImputeM %||% "nothing_to_impute",
      remove_imputed = isTRUE(input$remove_imputed),
      post_cor_filter = input$post_cor_filter,
      rsd_cutoff = filtered_corrected_r()$rsd_cutoff,
      remove_qc_average_pct_filter = isTRUE(input$remove_qc_average_pct_filter),
      transform = input$transform,
      ex_ISTD = isTRUE(input$ex_ISTD),
      keep_corrected_qcs = isTRUE(input$keep_corrected_qcs),
      corr_threshold = input$corr_threshold
    ))

    list(
      imputed            = imputed_r,
      corrected          = corrected_r,
      filtered_corrected = filtered_corrected_r,
      hotelling_res      = hotelling_res_r,
      transformed        = transformed_r,
      all_corr           = all_corr_r,
      params             = correct_params
    )
  })
}
