#----------- Reusable functions
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

#--------- 1.1 Upload Raw Data inputs
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

#---------- 1.2 Raw Data Inspection Inputs
#' Column selection for meta data
#' @keywords internal
#' @noRd
ui_nonmet_cols <- function(cols, ns = identity) {
  dropdown_choices <- c("Select a column..." = "", cols)
  
  tagList(
    htmltools::tags$h5("Select Non-Metabolite Columns"),
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

ui_control_class_selector <- function(df, ns) {
  classes <- unique(df$class[df$class != "QC"])
  dropdown_choices <- c("Select a class..." = "", classes)
  
  htmltools::tagList(
    htmltools::tags$hr(),
    htmltools::tags$h5("Select control class (optional)"),
    #htmltools::tags$p("If your data includes a control group, select it below. If not, check “No control group”."),
      bslib::tooltip(
        shiny::checkboxInput(ns("no_control"), "No control class", FALSE),
        "Check this if the dataset does not have a control class.",
        placement = "right"
      ),
      shiny::conditionalPanel(
        condition = sprintf("!input['%s']", ns("no_control")),
        bslib::tooltip(
          shiny::selectInput(
            ns("control_class"),
            "Control class",
            choices = dropdown_choices,
            selected = ""
          ),
          "Name of control samples in the class column. This class’s average is used to compute fold changes in the Excel file exported from this app. Fold changes are exported to a separate tab in the corrected-data Excel file.",
          placement = "right"
        )
      )
    )
}


#----------- 1.3 Filter Missing Values
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

#---------- 2.1 Choose Correction Settings inputs
#' Impute missing QC value options
#' @keywords internal
#' @noRd
ui_qc_impute <- function(df, metab_cols, ns = identity) {
  qc_df <- df %>% dplyr::filter(df$class == "QC")
  has_qc_na <- any(is.na(qc_df[, metab_cols]))
  
  if (has_qc_na) {
    
    label_with_info <- shiny::tagList(
      shiny::span("QC Sample Imputation Method"),
      bslib::popover(
        shiny::tags$button(
          type = "button",
          class = "btn btn-link p-0 ms-1",
          style = "text-decoration:none;",
          shiny::icon("circle-info")
        ),
        shiny::tagList(
          shiny::tags$p(
            "Choose how missing values in QC samples are imputed before correction."
          ),
          shiny::tags$ul(
            shiny::tags$li(shiny::strong("Metabolite median / mean: "),
                           "Across all samples."),
            shiny::tags$li(shiny::strong("QC-metabolite median / mean: "),
                           "Across QC samples only."),
            shiny::tags$li(shiny::strong("Minimum / half-minimum: "),
                           "Common for left-censored LC–MS data. Left-censored data occur when metabolite intensities fall below the instrument’s detection limit, so their exact values are unknown but known to be small."),
            shiny::tags$li(shiny::strong("KNN: "),
                           "k-nearest neighbors imputation."),
            shiny::tags$li(shiny::strong("Zero: "),
                           "Not recommended unless biologically justified.")
          )
        ),
        title = "QC imputation methods",
        placement = "right",
        options = list(
          container = "body",
          customClass = "popover-responsive"
        )
      )
    )
    
    shiny::radioButtons(
      ns("qcImputeM"),
      label = label_with_info,
      choices = list(
        "metabolite median" = "median",
        "metabolite mean" = "mean",
        "QC-metabolite median" = "class_median",
        "QC-metabolite mean" = "class_mean",
        "minimum value" = "min",
        "half minimum value" = "minHalf",
        "KNN" = "KNN",
        "zero" = "zero"
      ),
      selected = "median",
      inline = FALSE
    )
    
  } else {
    
    shiny::tags$div(
      shiny::icon("check-circle", class = "text-success"),
      shiny::span("No QC missing values")
    )
    
  }
}

#' Impute missing sample value options
#' @keywords internal
#' @noRd
ui_sample_impute <- function(df, metab_cols, ns = identity) {
  sam_df <- df %>% filter(df$class != "QC")
  has_sam_na <- any(is.na(sam_df[, metab_cols]))
  num_classes <- length(unique(sam_df$class))
  
  label_with_info <- shiny::tagList(
    shiny::span("Sample Imputation Method"),
    bslib::popover(
      shiny::tags$button(
        type = "button",
        class = "btn btn-link p-0 ms-1",
        style = "text-decoration:none;",
        shiny::icon("circle-info")
      ),
      shiny::tagList(
        shiny::tags$p(
          "Choose how missing values in samples are imputed before correction."
        ),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Metabolite median / mean: "),
                         "Across all samples."),
          shiny::tags$li(shiny::strong("Class-metabolite median / mean: "),
                         "Across samples grouping by class."),
          shiny::tags$li(shiny::strong("Minimum / half-minimum: "),
                         "Common for left-censored LC–MS data. Left-censored data occur when metabolite intensities fall below the instrument’s detection limit, so their exact values are unknown but known to be small."),
          shiny::tags$li(shiny::strong("KNN: "),
                         "k-nearest neighbors imputation."),
          shiny::tags$li(shiny::strong("Zero: "),
                         "Not recommended unless biologically justified.")
        )
      ),
      title = "Sample imputation methods",
      placement = "right",
      options = list(
        container = "body",
        customClass = "popover-responsive"
      )
    )
  )
  
  if (has_sam_na) {
    if (num_classes > 1) {
      radioButtons(
        inputId = ns("samImputeM"),
        label = label_with_info,
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
        label = label_with_info,
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

#' Correction method selection options
#' @keywords internal
#' @noRd
ui_correction_method <- function(df, ns = identity) {
  qc_per_batch <- df %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(qc_in_batch = sum(class == "QC", na.rm = TRUE), .groups = "drop")
  
  num_batches <- dplyr::n_distinct(df$batch)
  total_qcs   <- sum(df$class == "QC", na.rm = TRUE)
  
  label_with_info <- shiny::tagList(
    shiny::span("Correction Regression Model"),
    bslib::popover(
      shiny::tags$button(
        type = "button",
        class = "btn btn-link p-0 ms-1",
        style = "text-decoration:none;",
        shiny::icon("circle-info")
      ),
      report_text_correction_descriptions(),
      title = "What do these methods mean?",
      placement = "right",
      options = list(container = "body", customClass = "popover-responsive")
    )
  )
  
  # ---- Base choices depend ONLY on total number of QCs ----
  base_choices <- NULL
  base_selected <- NULL
  
  if (total_qcs <= 4) {
    base_choices <- list(
      "Local constant" = "LC",
      "Local linear"   = "LL"
    )
    base_selected <- "LL"
  } else if (total_qcs <= 8) {
    base_choices <- list(
      "Local constant"   = "LC",
      "Local linear"     = "LL",
      "Local polynomial" = "LOESS"
    )
    base_selected <- "LL"
  } else if (total_qcs <= 15) {
    base_choices <- list(
      "Local constant"   = "LC",
      "Local linear"     = "LL",
      "Local polynomial" = "LOESS",
      "Random forest"    = "RF"
    )
    base_selected <- "LOESS"
  } else {
    base_choices <- list(
      "Local constant"   = "LC",
      "Local linear"     = "LL",
      "Local polynomial" = "LOESS",
      "Random forest"    = "RF"
    )
    base_selected <- "RF"
  }
  
  # ---- Keep the existing batch-wise option logic unchanged ----
  choices <- base_choices
  selected <- base_selected
  
  #if (num_batches > 1 && !any(qc_per_batch$qc_in_batch < 5)) {
    #choices <- c(
      #choices,
      #list(
        #"Batchwise random forest"              = "BW_RF",
        #"Batchwise local polynomial fit (LOESS)" = "BW_LOESS"
      #)
    #)
    # keep selected as the base-selected method
  #}
  
  shiny::radioButtons(
    inputId  = ns("corMethod"),
    label    = label_with_info,
    choices  = choices,
    selected = selected
  )
}


#---------- 2.2 Post-Correction Filtering inputs
#' Post-correction filtering
#' @keywords internal
#' @noRd
ui_post_cor_filter <- function(ns) {
  shiny::tagList(
    tooltip(
      checkboxInput(ns("remove_imputed"), "Remove imputed values after correction", FALSE),
      "Check this box if you want to the corrected data to have the same missing values as the raw data.", 
      placement = "right"
    ),
    htmltools::tags$h5("QC RSD Filtering"),
    shiny::tags$div(
      style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
      shiny::tags$strong("RSD calculation"),
      bslib::popover(
        shiny::tags$button(
          type = "button",
          class = "btn btn-link p-0",
          style = "text-decoration:none;",
          shiny::icon("circle-info")
        ),
        report_text_rsd_cal(),
        title = "RSD% calculation",
        placement = "auto",
        options = list(container = "body",
                       customClass = "popover-responsive") 
      )
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
    ),
    htmltools::tags$h5("Candidate Extreme Values"),
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
}

#---------- 2.3 Post-Correction Transformation
#' Post-correction transformation
#' @keywords internal
#' @noRd
ui_post_cor_transform <- function(df, metab_cols, ns = identity) {
  has_istd <- any(grepl("ISTD|ITSD", metab_cols, ignore.case = FALSE))
  
  choices <- if (has_istd) {
    list(
      "Internal Standard Normalization" = "ISTD_norm",
      "Total Ratiometrically Normalized (TRN)" = "TRN",
      "None" = "none"
    )
  } else {
    list(
      "Total Ratiometric Normalization (TRN)" = "TRN",
      "None" = "none"
    )
  }
  
  label_with_info <- shiny::tagList(
    shiny::span("Method"),
    bslib::popover(
      shiny::tags$button(
        type = "button",
        class = "btn btn-link p-0 ms-1",
        style = "text-decoration:none;",
        shiny::icon("circle-info")
      ),
      report_text_transform_methods(),
      title = "Transformation methods",
      placement = "right",
      options = list(container = "body", customClass = "popover-responsive")
    )
  )
  
  shiny::tagList(
    shiny::radioButtons(
      ns("transform"),
      label = label_with_info,
      choices = choices,
      selected = "none"
    ),
    bslib::tooltip(
      shiny::checkboxInput(
        ns("ex_ISTD"),
        "Exclude internal standards from post-correction transformation.",
        TRUE
      ),
      "Check this box if you do not want internal standards to be included in the transformation calculation. Internal standards will appear in this table, but not in the '3. Scaled or Normalized' tab of the Excel file.",
      placement = "right"
    ),
    shiny::conditionalPanel(
      condition = sprintf("input['%s'] === 'TRN'", ns("transform")),
      bslib::tooltip(
        shiny::checkboxInput(ns("trn_withhold_checkbox"), "Withhold column(s) from TRN", FALSE),
        "Check this box if there are any columns that should not count in TRN (i.e. TIC column). Sample, batch, class and order are already excluded.",
        placement = "right"
      )
    )
  )
}

#----------- 2.4 Metabolite Correlations inputs
#' metabolite correlation slider
#' @keywords internal
#' @noRd
ui_correlation_slider <- function(ns) {
  tooltip(
    sliderInput(ns("corr_threshold"), "Pearson's r range", 0.9, 1, value = c(0.99, 1), step = 0.005),
    "Pairs of metabolites with Pearson's r within this range will be displayed on the rigth after clicking the 'Compute Metabolite Correlations' button.", 
    placement = "right"
  )
}

#---------- 3.1 Scatter Plot Evaluation input

#---------- 3.2 RSD Evaluation input
#' visualization rsd evaluation
#' @keywords internal
#' @noRd
ui_rsd_eval <- function(ns) {
  tagList(
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

#---------- 3.3 PCA Evaluation input
#' visualization pca evaluation
#' @keywords internal
#' @noRd
ui_pca_eval <- function(ns){
  tagList(
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

#----------- 3.4 Select Figure Format input
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