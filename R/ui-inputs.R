#' File upload
#' @keywords internal
#' @noRd
ui_file_upload <- function(ns) {
    fileInput(
      ns("file1"),
      "Choose Raw Data File (.csv, .xls, or .xlsx)",
      accept = c(".csv", ".xls", ".xlsx"),
      buttonLabel = "Browse...",
      placeholder = "No file selected"
    )
}

#' Reusable titled sidebar
#' @keywords internal
#' @noRd
ui_sidebar_block <- function(title, ..., help = NULL, width = 400, position = "left") {
  bslib::sidebar(
    tags$h4(title),
    ...,
    if (!is.null(help)) lapply(help, tags$h6),
    width = width,
    position = position
  )
}

#' Scrollable table wrapper
#' @keywords internal
#' @noRd
ui_table_scroll <- function(outputId, ns, height = "400px") {
  div(
    style = paste0("overflow:auto; max-height:", height, ";"),
    tableOutput(ns(outputId))
  )
}

#' Column selection for meta data
#' @keywords internal
#' @noRd
ui_nonmet_cols <- function(cols, ns = identity) {
  dropdown_choices <- c("Select a column..." = "", cols)
  
  tagList(
    tooltip(
      selectInput(ns("sample_col"), "sample column", dropdown_choices, ""),
      "Column that contains unique sample names.",
      placement ="right"
    ),
    tooltip(
      checkboxInput(ns("single_batch"), "no batch column", FALSE), 
      "check this box if your raw data does not have a column indicating batch. All samples will be assigned the same batch for correction.",
      placement = "right"
    ),
    conditionalPanel(
      condition = sprintf("!input['%s']", ns("single_batch")),
      tooltip(
        selectInput(ns("batch_col"), "batch column", dropdown_choices, ""),
        "Column that contains batch information.",
        placement ="right"
      )
    ),
    tooltip(
      selectInput(ns("class_col"), "class column", dropdown_choices, ""),
      "Column that indicates the type of sample. Must contain QC samples labeled as 'NA', 'QC', 'Qc', or 'qc'.",
      placement ="right"
    ),
    tooltip(
      selectInput(
        ns("order_col"),
        "injection order column",
        dropdown_choices,
        ""
      ),
      "Column that indicates injection order.",
      placement = "right"
    )
  )
}

#' Toggle for withholding extra columns
#' @keywords internal
#' @noRd
ui_withhold_toggle <- function(ns) {
  tooltip(
    checkboxInput(ns("withhold_cols"),
                  "Withhold additional columns from correction", FALSE),
    "Select if there are extra non-metabolite or specific metabolite columns to withhold.",
    placement = "right"
  )
}

#' Count input for how many columns to withhold
#' @keywords internal
#' @noRd
ui_withhold_count <- function(ns, max_withhold) {
  numericInput(ns("n_withhold"),
               "Number of columns to withhold",
               value = if (max_withhold > 0) 1 else 0,
               min   = 0,
               max   = max_withhold,
               step  = 1)
}

#' metabolite correlation slider
#' @keywords internal
#' @noRd
ui_corr_slider <- function(ns) {
  tooltip(
    sliderInput(ns("corr_threshold"), "Pearson's r range", 0.9, 1, value = c(0.99, 1), step = 0.005),
    "Pairs of metabolites with Pearson's R within this range will be displayed here.", 
    placement = "right"
  )
}
#' metabolite correlation slider in corrected or transformed data
#' @keywords internal
#' @noRd
ui_tc_corr_slider <- function(ns) {
  tagList(tooltip(
    radioButtons(
      ns("tc_corr_data"),
      "Compute metabolite correlations for",
      list(
        "Corrected data" = "filtered_cor_data",
        "Transformed and corrected data" = "transformed_cor_data"
      ),
      "filtered_cor_data"
    ),
    "all pairwise metabolite correlations will be computed in the data set you select.",
    placement = "right"
  ),
  tooltip(
    sliderInput(ns("tc_corr_threshold"), "Pearson's r range", 0.9, 1, value = c(0.99, 1), step = 0.005),
    "Pairs of metabolites with Pearson's R within this range will be displayed here.", 
    placement = "right"
  )
  )
}

#' missing value filter slider
#' @keywords internal
#' @noRd
ui_filter_slider <- function(ns) {
  tooltip(
    sliderInput(ns("mv_cutoff"), "Acceptable % missing per metabolite", 0, 100, 20),
    "Metabolites with missing % above this threshold are removed.", 
    placement = "right"
  )
}

#' Repeated selectors for which columns to withhold
#' @param ids character vector of input ids to render (e.g., "withhold_col_1")
#' @param cols candidate column names
#' @param prev named character of previous selections for each id (same length as ids)
#' @keywords internal
#' @noRd
ui_withhold_selectors <- function(ids, cols, prev, ns) {
  if (!length(ids)) return(NULL)
  # keep uniqueness across the repeated selects
  lapply(seq_along(ids), function(i) {
    id    <- ids[i]
    prior <- prev[[i]] %||% ""
    other <- setdiff(prev, prior)
    choices_i <- c("Select a column..." = "", setdiff(cols, other))
    selectInput(
      ns(id),
      label   = paste("Select column to withhold #", i),
      choices = choices_i,
      selected = if (nzchar(prior) && prior %in% choices_i) prior else ""
    )
  })
}

#' Impute missing QC value options for section 2.1 Choose Correction settings
#' @keywords internal
#' @noRd
ui_qc_impute <- function(df, metab_cols, ns = identity) {
  qc_df <- df %>% filter(df$class == "QC")
  has_qc_na <- any(is.na(qc_df[, metab_cols]))
  
  if (has_qc_na) {
    radioButtons(
      ns("qcImputeM"),
      "QC Sample Imputation Method",
      list(
        "metabolite median" = "median",
        "metabolite mean" = "mean",
        "QC-metabolite median" = "class_median",
        "QC-metabolite mean" = "class_mean",
        "minimum value" = "min",
        "half minimum value" = "minHalf",
        "KNN" = "KNN",
        "zero" = "zero"
      ),
      "median",
      FALSE
    )
  } else {
    tags$div(icon("check-circle", class = "text-success"),
             span("No QC missing values"))
  }
}

#' Impute missing sample value options for section 2.1 Choose Correction settings
#' @keywords internal
#' @noRd
ui_sample_impute <- function(df, metab_cols, ns = identity) {
  sam_df <- df %>% filter(df$class != "QC")
  has_sam_na <- any(is.na(sam_df[, metab_cols]))
  num_classes <- length(unique(sam_df$class))
  
  if (has_sam_na) {
    if (num_classes > 1) {
      radioButtons(
        inputId = ns("samImputeM"),
        label = "Sample Imputation Method",
        choices = list(
          "metabolite median" = "median",
          "metabolite mean" = "mean",
          "class-metabolite median" = "class_median",
          "class-metabolite mean" = "class_mean",
          "minimum value" = "min",
          "half minimum value" = "minHalf",
          "KNN" = "KNN",
          "zero" = "zero"
        ),
        selected = "median",
        inline = FALSE
      )
    } else {
      radioButtons(
        inputId = ns("samImputeM"),
        label = "Sample Imputation Method",
        choices = list(
          "metabolite median" = "median",
          "metabolite mean" = "mean",
          "minimum value" = "min",
          "half minimum value" = "minHalf",
          "KNN" = "KNN",
          "zero" = "zero"
        ),
        selected = "median",
        inline = FALSE
      )
    }
  } else {
    tags$div(icon("check-circle", class = "text-success"),
             span("No Sample missing values"))
  }
}

#' Correction method selection options for section 2.1 Choose Correction settings
#' @keywords internal
#' @noRd
ui_correction_method <- function(df, ns = identity) {
  qc_per_batch <- df %>%
    group_by(batch) %>%
    summarise(qc_in_batch = sum(class == "QC"), .groups = "drop")
  num_batches <- length(unique(df$batch))
  
  if (num_batches == 1) {
    if (any(qc_per_batch$qc_in_batch <= 5)) {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list("Local Polynomial Fit (LOESS)" = "LOESS"),
        selected = "LOESS"
      )
    } else {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list(
          "Random Forest" = "RF",
          "Local Polynomial Fit (LOESS)" = "LOESS"
        ),
        selected = "RF"
      )
    }
  } else {
    if (any(qc_per_batch$qc_in_batch < 5)) {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list(
          "Random Forest" = "RF",
          "Local Polynomial Fit (LOESS)" = "LOESS"
        ),
        selected = "RF"
      )
    } else {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list(
          "Random Forest" = "RF",
          "Local Polynomial Fit (LOESS)" = "LOESS",
          "Batchwise Random Forest" = "BW_RF",
          "Batchwise Local polynomial fit (LOESS)" = "BW_LOESS"
        ),
        selected = "RF"
      )
    }
  }
}

#' Post-correction filtering
#' @keywords internal
#' @noRd
ui_post_cor_filter <- function(ns) {
  tagList(
    tooltip(
      checkboxInput(ns("remove_imputed"), "Remove imputed values after correction", FALSE),
      "Check this box if you want to the corrected data to have the same missing values as the raw data.", 
      placement = "right"
    ),
    tooltip(
      checkboxInput(ns("post_cor_filter"), "Don't filter metabolites based on QC RSD%", FALSE),
      "Check this box if you don't want any metabolites removed post-correction.", 
      placement = "right"
    ),
    conditionalPanel(
      condition = sprintf("!input['%s']", ns("post_cor_filter")),
      tooltip(
       sliderInput(ns("rsd_filter"),"Metabolite RSD% threshold for QC samples", 0, 100, 20),
        "Metabolites with QC RSD% above this value will be removed from the corrected data.", 
       placement = "right"
      )
    )
  )
}

#' Post-correction transformation
#' @keywords internal
#' @noRd
ui_post_cor_transform <- function(df, metab_cols, ns = identity) {
  has_istd <- any(grepl("^(ISTD|ITSD)", metab_cols, ignore.case = TRUE))
  
  choices <- if (has_istd) {
    list(
      "Internal Standard Normalization" = "ISTD_norm",
      "Total Ratiometically Normalized (TRN)" = "TRN",
      "None" = "none"
    )
  } else {
    list(
      "Total Ratiometic Normalization (TRN)" = "TRN",
      "None" = "none"
    )
  }
  
  tagList(
    radioButtons(
      ns("transform"),
      "Method",
      choices  = choices,
      selected = "none"
    ),
    tooltip(
      checkboxInput(
        ns("ex_ISTD"),
        "Exclude Internal Standards from post-correction transformation/normalization.",
        TRUE
      ),
      "Check this box if you do not want internal standards to be included in the transformation or normalization calculation.",
      placement = "right"
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] === 'TRN'", ns("transform")),
      tooltip(
        checkboxInput(ns("trn_withhold_checkbox"), "Withhold column(s) from TRN", FALSE),
        "Check this box if there are any columns that should not count in TRN (i.e. TIC column). Sample, batch, class and order are already excluded.",
        placement = "right"
      )
    )
  )
}

#' Options for outlier detection
#' @keywords internal
#' @noRd
ui_detect_outliers_options <- function(ns) {
  tooltip(
    radioButtons(
      ns("out_data"),
      "Detect extreme values in",
      list(
        "Corrected data" = "filtered_cor_data",
        "Transformed and corrected data" = "transformed_cor_data"
      ),
      "filtered_cor_data"
    ),
    "Potential extreme values will be detected in the data set you select.",
    placement = "right"
  )
}

#' visualization rsd evaluation
#' @keywords internal
#' @noRd
ui_rsd_eval <- function(ns) {
  tagList(
    tags$h6("Evaluate correction method by the change in relative standard deviation (RSD)."),
    radioButtons(ns("rsd_plot_type"),
                 "Visualize Changes in RSD by",
                 list("Distribution" = "dist",
                      "Scatter Plot" = "scatter")
    ),
    radioButtons(ns("rsd_compare"), 
                 "Compare raw data to", 
                 list("Corrected data" = "filtered_cor_data", 
                      "Transformed and corrected data" = "transformed_cor_data"), 
                 "filtered_cor_data"),
    radioButtons(ns("rsd_cal"), 
                 "Calculate RSD by", 
                 list("Metabolite" = "met", "Class and Metabolite" = "class_met"),
                 "met")
  )
}

#' visualization pca evaluation
#' @keywords internal
#' @noRd
ui_pca_eval <- function(ns){
  tagList(
    tags$h6("Evaluate correction using principal component analysis (PCA)."),
    radioButtons(ns("pca_compare"), 
                 "Compare raw data to", 
                 list("Corrected data" = "filtered_cor_data", 
                      "Transformed and corrected data" = "transformed_cor_data"), 
                 "filtered_cor_data"),
    radioButtons(ns("color_col"), 
                 "Color PCA by", 
                 list("batch" = "batch", "class" = "class", "order" = "order"), 
                 "batch")
  )
}

#' Visualization downloading figure format
#' @keywords internal
#' @noRd
ui_fig_format <- function(ns) {
  tooltip(
    radioButtons(ns("fig_format"), "Select figure format:", c("PDF" = "pdf", "PNG" = "png"), "pdf"),
    "All figures will be saved in this format after clicking download button here or on tab 4. Export Corrected Data, Plots, and Report", 
    placement = "right"
  )
}