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
    htmltools::tags$h5("Select Required Metadata Columns"),
    tooltip(
      selectInput(ns("sample_col"), "sample column", dropdown_choices, ""),
      "Column that sample names. Data cannot have repeated sample names.",
      placement = "right"
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
        placement = "right"
      )
    ),
    tooltip(
      selectInput(ns("class_col"), "class column", dropdown_choices, ""),
      "Column that indicates the type of sample. Must contain QC samples labeled as 'NA', 'QC', 'Qc', or 'qc'.",
      placement = "right"
    ),
    tooltip(
      selectInput(
        ns("order_col"),
        "injection order column",
        dropdown_choices,
        ""
      ),
      "Column that indicates the order in which samples were injected into the instrument.",
      placement = "right"
    )
  )
}

#' Toggle for withholding extra columns
#' @keywords internal
#' @noRd
ui_withhold_toggle <- function(ns) {
  tooltip(
    checkboxInput(
      ns("withhold_cols"),
      "Data contains additional metadata columns", FALSE
    ),
    "Select if there are extra non-metabolite columns in the dataset.",
    placement = "right"
  )
}

ui_control_class_selector <- function(df, ns) {
  classes <- unique(df$class[df$class != "QC"])
  dropdown_choices <- c("Select a class..." = "", classes)

  htmltools::tagList(
    htmltools::tags$hr(),
    htmltools::tags$h5("Select control class (optional)"),
    # htmltools::tags$p("If your data includes a control group, select it below. If not, check “No control group”."),
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


#----------- 1.3 Raw Data Filtering
#' slider for blank threshold and checkbox to optionally filter metabolites
#' that fail the threshold.
#' @keywords internal
#' @noRd
ui_blank_threshold_controls <- function(ns = identity,
                                        threshold = 3,
                                        remove_default = FALSE) {
  shiny::tagList(
    shiny::tags$div(
      style = "display:flex; align-items:center; justify-content:space-between; gap: 8px; margin-bottom: 8px;",
      shiny::tags$strong("Blank Threshold Filter"),
      bslib::popover(
        shiny::tags$button(
          type = "button",
          class = "btn btn-link p-0",
          style = "text-decoration:none;",
          shiny::icon("circle-info")
        ),
        shiny::tags$div(
          shiny::tags$p(
            "For each metabolite, the average QC intensity is compared with the average blank or processing blank intensity."
          ),
          shiny::tags$p(
            "A metabolite fails when its QC average is less than the selected multiplier times the blank average."
          ),
          shiny::tags$p(
            "Failed metabolites can be flagged only, or removed before missing-value filtering."
          )
        ),
        title = "Why filter metabolites based on blanks?",
        placement = "auto",
        options = list(
          container = "body",
          customClass = "popover-responsive"
        )
      )
    ),
    shiny::sliderInput(
      inputId = ns("blank_threshold"),
      label = "Blank threshold multiplier",
      min = 1,
      max = 20,
      value = threshold,
      step = 1
    ),
    shiny::checkboxInput(
      inputId = ns("remove_blank_threshold_cols"),
      label = "Remove metabolites that fail blank threshold",
      value = remove_default
    ),
    shiny::tags$hr()
  )
}

#' missing value filter slider
#' @keywords internal
#' @noRd
ui_filter_slider <- function(ns) {
  tooltip(
    sliderInput(ns("mv_cutoff"), "Acceptable % missing per metabolite", 0, 100, 20),
    "Metabolites with missing % above this threshold for at least 1 class are removed.",
    placement = "right"
  )
}

#---------- 2.1 Choose Correction Settings inputs
#' Impute missing QC value options
#' @keywords internal
#' @noRd
ui_qc_impute <- function(df, metab_cols, ns = identity) {
  qc_df <- df |>
    dplyr::filter(df$class == "QC")
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
            shiny::tags$li(
              shiny::strong("Metabolite median / mean: "),
              "Across all samples."
            ),
            shiny::tags$li(
              shiny::strong("QC-metabolite median / mean: "),
              "Across QC samples only."
            ),
            shiny::tags$li(
              shiny::strong("Minimum / half-minimum: "),
              "Common for left-censored LC–MS data. Left-censored data occur when metabolite intensities fall below the instrument’s detection limit, so their exact values are unknown but known to be small."
            ),
            shiny::tags$li(
              shiny::strong("KNN: "),
              "k-nearest neighbors imputation."
            ),
            shiny::tags$li(
              shiny::strong("Zero: "),
              "Not recommended unless biologically justified."
            )
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
  sam_df <- df |>
    dplyr::filter(df$class != "QC")
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
          shiny::tags$li(
            shiny::strong("Metabolite median / mean: "),
            "Across all samples."
          ),
          shiny::tags$li(
            shiny::strong("Class-metabolite median / mean: "),
            "Across samples grouping by class."
          ),
          shiny::tags$li(
            shiny::strong("Minimum / half-minimum: "),
            "Common for left-censored LC–MS data. Left-censored data occur when metabolite intensities fall below the instrument’s detection limit, so their exact values are unknown but known to be small."
          ),
          shiny::tags$li(
            shiny::strong("KNN: "),
            "k-nearest neighbors imputation."
          ),
          shiny::tags$li(
            shiny::strong("Zero: "),
            "Not recommended unless biologically justified."
          )
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
    tags$div(
      icon("check-circle", class = "text-success"),
      span("No Sample missing values")
    )
  }
}

#' Correction method selection options
#' @keywords internal
#' @noRd
ui_correction_method <- function(df, ns = identity) {
  stopifnot(is.data.frame(df))

  required_cols <- c("class", "batch", "order")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "df is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  qc_per_batch <- df |>
    dplyr::group_by(batch) |>
    dplyr::summarise(
      qc_in_batch = sum(class == "QC", na.rm = TRUE),
      .groups = "drop"
    )

  num_batches <- dplyr::n_distinct(df$batch)
  total_qcs <- sum(df$class == "QC", na.rm = TRUE)

  # Compute QC spacing from injection order
  qc_gap_stats <- NULL
  ord <- df$order
  is_qc <- df$class == "QC"

  keep <- !is.na(ord) & !is.na(is_qc)
  ord <- ord[keep]
  is_qc <- is_qc[keep]

  if (!is.numeric(ord)) {
    ord_num <- suppressWarnings(as.numeric(ord))

    if (!all(is.na(ord_num))) {
      ord <- ord_num
    }
  }

  if (is.numeric(ord) && sum(is_qc, na.rm = TRUE) >= 2L) {
    qc_orders <- sort(ord[is_qc])
    gaps <- diff(qc_orders)

    if (length(gaps) > 0L) {
      qc_gap_stats <- list(
        min_gap = min(gaps, na.rm = TRUE),
        median_gap = stats::median(gaps, na.rm = TRUE),
        max_gap = max(gaps, na.rm = TRUE)
      )
    }
  }

  median_gap <- qc_gap_stats$median_gap %||% NA_real_
  max_gap <- qc_gap_stats$max_gap %||% NA_real_

  has_spacing <- !is.null(qc_gap_stats) &&
    is.finite(median_gap) &&
    is.finite(max_gap)

  # QC-count eligibility
  allow_lc <- total_qcs >= 1L
  allow_ll <- total_qcs >= 3L
  allow_loess <- total_qcs >= 5L
  allow_rf <- total_qcs >= 9L

  # Spacing support: used for recommendation only
  supports_loess <- has_spacing &&
    median_gap <= 10 &&
    max_gap <= 15

  supports_rf <- has_spacing &&
    median_gap <= 9 &&
    max_gap <= 10

  # All methods remain visible.
  # `allowed` controls whether the radio option is enabled.
  method_options <- tibble::tibble(
    label = c(
      "Local constant",
      "Local linear",
      "Local polynomial",
      "Random forest"
    ),
    value = c("LC", "LL", "LOESS", "RF"),
    allowed = c(
      allow_lc,
      allow_ll,
      allow_loess,
      allow_rf
    ),
    unavailable_reason = c(
      "Requires at least 1 QC sample.",
      "Requires at least 3 QC samples.",
      "Requires at least 5 QC samples.",
      "Requires at least 9 QC samples."
    )
  )

  # Recommended selection
  selected <- if (allow_rf && supports_rf) {
    "RF"
  } else if (allow_loess && supports_loess) {
    "LOESS"
  } else if (allow_ll) {
    "LL"
  } else if (allow_lc) {
    "LC"
  } else {
    character(0)
  }

  available_values <- method_options$value[method_options$allowed]

  if (length(selected) == 0L || !selected %in% available_values) {
    selected <- if (length(available_values) > 0L) {
      available_values[[1L]]
    } else {
      character(0)
    }
  }

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
      options = list(
        container = "body",
        customClass = "popover-responsive"
      )
    )
  )

  input_id <- ns("corMethod")
  label_id <- paste0(input_id, "-label")

  make_radio_choice <- function(label, value, allowed, unavailable_reason) {
    is_selected <- length(selected) == 1L && identical(value, selected)

    shiny::tags$div(
      class = if (isTRUE(allowed)) "radio" else "radio disabled",
      title = if (isTRUE(allowed)) NULL else unavailable_reason,
      shiny::tags$label(
        shiny::tags$input(
          type = "radio",
          name = input_id,
          value = value,
          checked = if (is_selected) "checked" else NULL,
          disabled = if (!isTRUE(allowed)) "disabled" else NULL
        ),
        shiny::tags$span(label)
      )
    )
  }

  radio_choices <- purrr::pmap(
    list(
      label = method_options$label,
      value = method_options$value,
      allowed = method_options$allowed,
      unavailable_reason = method_options$unavailable_reason
    ),
    make_radio_choice
  )

  shiny::tagList(
    shiny::tags$style(
      htmltools::HTML(
        "
        .cor-method-radio .radio.disabled label {
          color: #6c757d;
          opacity: 0.65;
          cursor: not-allowed;
        }

        .cor-method-radio .radio.disabled input {
          cursor: not-allowed;
        }
        "
      )
    ),
    shiny::tags$div(
      id = input_id,
      class = "form-group shiny-input-radiogroup shiny-input-container cor-method-radio",
      role = "radiogroup",
      `aria-labelledby` = label_id,
      shiny::tags$label(
        class = "control-label",
        id = label_id,
        `for` = input_id,
        label_with_info
      ),
      shiny::tags$div(
        class = "shiny-options-group",
        radio_choices
      )
    )
  )
}


#---------- 2.2 Post-Correction Filtering inputs
#' Post-correction filtering
#' @keywords internal
#' @noRd
ui_post_cor_filter <- function(ns) {
  shiny::tagList(
    htmltools::tags$h5("Imputed Values"),
    tooltip(
      shiny::checkboxInput(
        inputId = ns("remove_imputed"),
        label = "Remove imputed values after correction",
        value = FALSE
      ),
      "Check this box if you want the corrected data to have the same missing values as the raw data.",
      placement = "right"
    ),
    shiny::tags$hr(),
    htmltools::tags$h5("Sample Distance From QC Filtering"),
    tooltip(
      shiny::sliderInput(
        inputId = ns("qc_average_pct_threshold"),
        label = "Minimum % distance from QC average",
        min = 0,
        max = 200,
        value = 50,
        step = 5,
        post = "%"
      ),
      "Metabolites are flagged when the average non-QC sample intensity differs from the average QC intensity by at least the selected percentage threshold.",
      placement = "right"
    ),
    tooltip(
      shiny::checkboxInput(
        inputId = ns("remove_qc_average_pct_filter"),
        label = "Remove metabolites far from QC average",
        value = FALSE
      ),
      "Removes metabolites where the average non-QC sample intensity differs from the average QC intensity by at least the selected percentage.",
      placement = "right"
    ),
    shiny::tags$hr(),
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
        options = list(
          container = "body",
          customClass = "popover-responsive"
        )
      )
    ),
    tooltip(
      shiny::checkboxInput(
        inputId = ns("post_cor_filter"),
        label = "Don't filter metabolites based on QC RSD%",
        value = FALSE
      ),
      "Check this box if you don't want any metabolites removed post-correction.",
      placement = "right"
    ),
    shiny::conditionalPanel(
      condition = sprintf("!input['%s']", ns("post_cor_filter")),
      tooltip(
        shiny::sliderInput(
          inputId = ns("rsd_filter"),
          label = "Metabolite RSD% threshold for QC samples",
          min = 0,
          max = 100,
          value = 30,
          step = 5
        ),
        "Metabolites with QC RSD% above this value will be removed from the corrected data.",
        placement = "right"
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
      "Probabilistic Quotient Normalization (PQN)" = "PQN",
      "None" = "none"
    )
  } else {
    list(
      "Total Ratiometric Normalization (TRN)" = "TRN",
      "Probabilistic Quotient Normalization (PQN)" = "PQN",
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
    radioButtons(
      ns("rsd_plot_type"),
      "Visualize Changes in RSD by",
      list(
        "Distribution" = "dist",
        "Scatter Plot" = "scatter"
      )
    ),
    radioButtons(
      ns("rsd_compare"),
      "Compare raw data to",
      list(
        "Corrected data" = "filtered_cor_data",
        "Transformed and corrected data" = "transformed_cor_data"
      ),
      "filtered_cor_data"
    ),
    radioButtons(
      ns("rsd_cal"),
      "Calculate RSD by",
      list("Metabolite" = "met", "Class and Metabolite" = "class_met"),
      "met"
    )
  )
}

#---------- 3.3 PCA Evaluation input
#' visualization pca evaluation
#' @keywords internal
#' @noRd
ui_pca_eval <- function(meta_df, ns) {
  color_options <- setdiff(names(meta_df), "sample")
  shape_options <- setdiff(names(meta_df), c("sample", "order"))

  # optional prioritization
  priority <- c("batch", "class", "order")
  color_options <- c(
    intersect(priority, color_options),
    setdiff(color_options, priority)
  )
  shape_options <- c(
    intersect(c("batch", "class"), shape_options),
    setdiff(shape_options, c("batch", "class"))
  )

  tagList(
    radioButtons(
      ns("pca_compare"),
      "Compare raw data to",
      choices = c(
        "Corrected data" = "filtered_cor_data",
        "Transformed and corrected data" = "transformed_cor_data"
      ),
      selected = "filtered_cor_data"
    ),
    radioButtons(
      ns("color_col"),
      "Color PCA by",
      choices = stats::setNames(color_options, color_options),
      selected = color_options[1]
    ),
    radioButtons(
      ns("shape_col"),
      "Sample shape defined by",
      choices = stats::setNames(shape_options, shape_options),
      selected = shape_options[1]
    )
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
