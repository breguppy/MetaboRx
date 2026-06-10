#' @keywords internal
#' @noRd

# Column-warning text for section 1.2 Select non-metabolite columns
ui_column_warning <- function(data, selected) {
  warnings <- list()

  if (any(selected == "")) {
    warnings[[length(warnings) + 1]] <- tags$span(
      style = "color:darkorange; font-weight:bold;",
      icon("exclamation-triangle"),
      " Please select all four columns."
    )
  } else if (length(unique(selected)) < 4) {
    warnings[[length(warnings) + 1]] <- tags$span(
      style = "color:darkred; font-weight:bold;",
      icon("exclamation-triangle"),
      " The same column is selected multiple times."
    )
  }

  if (!any(selected == "") && length(unique(selected)) == 4) {
    samp_vec <- data[[selected[1]]]
    order_vec <- data[[selected[4]]]

    if (anyDuplicated(samp_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(
        style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
        icon("exclamation-triangle"),
        " Duplicate sample names detected!"
      )
    }
    if (!is.numeric(order_vec)) {
      warnings[[length(warnings) + 1]] <- tags$span(
        style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
        icon("exclamation-triangle"),
        " Order column must contain numbers."
      )
    }
    if (anyDuplicated(order_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(
        style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
        icon("exclamation-triangle"),
        " Duplicate order values detected!"
      )
    }
  }

  if (length(warnings) == 0) {
    return(NULL)
  } else if (length(warnings) == 1) {
    return(warnings[[1]])
  } else {
    return(do.call(tagList, warnings))
  }
}

# QC missing value warning for section 2.1 Choose Correction settings
ui_qc_missing_warning <- function(df) {
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  qc_idx <- which(df$class == "QC")
  n_missv <- sum(is.na(df[qc_idx, metab_cols]))

  if (n_missv > 0) {
    tags$span(
      style = "color:darkred; font-weight:bold;",
      icon("exclamation-triangle"),
      paste(" ", n_missv, " values missing from QC samples")
    )
  } else {
    NULL
  }
}


ui_how_to_correct <- function(df,
                              qc_label = "QC",
                              class_col = "class",
                              order_col = "order") {
  stopifnot(is.data.frame(df))
  stopifnot(class_col %in% names(df))

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  total_qcs <- sum(df[[class_col]] == qc_label, na.rm = TRUE)

  # Compute QC spacing using injection order
  qc_gap_stats <- NULL
  if (order_col %in% names(df)) {
    ord <- df[[order_col]]
    is_qc <- df[[class_col]] == qc_label

    keep <- !is.na(ord) & !is.na(is_qc)
    ord <- ord[keep]
    is_qc <- is_qc[keep]

    if (!is.numeric(ord)) {
      ord_num <- suppressWarnings(as.numeric(ord))
      if (!all(is.na(ord_num))) {
        ord <- ord_num
      }
    }

    if (is.numeric(ord) && sum(is_qc) >= 2L) {
      qc_orders <- sort(ord[is_qc])
      gaps <- diff(qc_orders)
      qc_gap_stats <- list(
        min_gap = min(gaps, na.rm = TRUE),
        median_gap = stats::median(gaps, na.rm = TRUE),
        max_gap = max(gaps, na.rm = TRUE)
      )
    }
  }

  gap_line <- if (!is.null(qc_gap_stats)) {
    sprintf(
      "QC spacing (injection order): median gap = %s, max gap = %s",
      format(qc_gap_stats$median_gap, digits = 3),
      format(qc_gap_stats$max_gap, digits = 3)
    )
  } else {
    "QC spacing (injection order): unavailable (need \u22652 QCs with valid order)"
  }

  summary_bits <- htmltools::tagList(
    htmltools::tags$div(
      class = "small text-muted",
      sprintf("Total QCs detected: %d", total_qcs)
    ),
    htmltools::tags$div(
      class = "small text-muted",
      gap_line
    )
  )

  warn_span <- function(show, text) {
    if (!isTRUE(show)) {
      return(NULL)
    }

    htmltools::tags$span(
      shiny::icon("exclamation-triangle", class = "text-danger-emphasis"),
      htmltools::tags$span(style = "margin-left: 6px;", text),
      style = "display:block; margin-top:4px;"
    )
  }

  has_spacing <- !is.null(qc_gap_stats)
  median_gap <- qc_gap_stats$median_gap %||% NA_real_
  max_gap <- qc_gap_stats$max_gap %||% NA_real_

  # Spacing categories
  # Random forest gets its own stricter threshold:
  # median gap must be <= 9, and max gap must be reasonably controlled.
  spacing_dense_rf <- has_spacing &&
    is.finite(median_gap) && is.finite(max_gap) &&
    median_gap <= 9 && max_gap <= 10

  spacing_reasonable_poly <- has_spacing &&
    is.finite(median_gap) && is.finite(max_gap) &&
    median_gap <= 10 && max_gap <= 15

  spacing_sparse <- has_spacing &&
    is.finite(max_gap) &&
    max_gap > 15

  # Method eligibility based on QC count only
  allow_local_constant <- total_qcs >= 1
  allow_local_linear <- total_qcs >= 3
  allow_local_polynomial <- total_qcs >= 5
  allow_random_forest <- total_qcs >= 9

  recommendation <- NULL
  rationale <- NULL

  if (total_qcs < 3) {
    recommendation <- "Local constant"
    rationale <- paste(
      "Very limited QC support. Use the least flexible method.",
      "Higher-flexibility methods should not be used at this QC count."
    )
  } else if (total_qcs <= 4) {
    if (spacing_sparse) {
      recommendation <- "Local constant"
      rationale <- paste(
        "Few QCs and wide spacing between them.",
        "Use the most conservative method."
      )
    } else {
      recommendation <- "Local linear"
      rationale <- paste(
        "Few QCs, but enough support for a modestly flexible fit.",
        "Local linear is the highest-supported method in this range."
      )
    }
  } else if (total_qcs <= 8) {
    if (spacing_reasonable_poly) {
      recommendation <- "Local polynomial (QC-RLSC)"
      rationale <- paste(
        "Moderate QC count with acceptable spacing.",
        "A moderately flexible smoother is supported."
      )
    } else {
      recommendation <- "Local linear"
      rationale <- paste(
        "Moderate QC count, but spacing is not tight enough for local polynomial.",
        "Prefer the more stable option."
      )
    }
  } else {
    if (spacing_dense_rf) {
      recommendation <- "Random forest (QC-RFSC)"
      rationale <- paste(
        "QC count is high enough and median QC spacing is tight enough to support the highest-flexibility method."
      )
    } else if (spacing_reasonable_poly) {
      recommendation <- "Local polynomial (QC-RLSC)"
      rationale <- paste(
        "QC support is good, but not strong enough for random forest.",
        "Use the intermediate-flexibility method."
      )
    } else {
      recommendation <- "Local linear"
      rationale <- paste(
        "Although QC count is adequate, spacing is too wide for more flexible methods.",
        "Prefer the more stable choice."
      )
    }
  }

  rec_block <- htmltools::tags$div(
    style = "margin-top: 10px; margin-bottom: 10px;",
    htmltools::tags$strong("Recommended method: "),
    recommendation,
    htmltools::tags$div(
      class = "small text-muted",
      style = "margin-top: 4px;",
      rationale
    )
  )

  items <- list()

  if (allow_local_constant) {
    items <- c(items, list(
      htmltools::tags$li(
        htmltools::tags$strong("Local constant: "),
        "Best when QC support is very limited or QCs are widely spaced. Least flexible, most conservative option."
      )
    ))
  }

  if (allow_local_linear) {
    items <- c(items, list(
      htmltools::tags$li(
        htmltools::tags$strong("Local linear: "),
        "Best when there are a small to moderate number of QCs. More flexible than local constant, but still relatively stable."
      )
    ))
  }

  if (allow_local_polynomial) {
    items <- c(items, list(
      htmltools::tags$li(
        htmltools::tags$strong("Local polynomial (QC-RLSC): "),
        "Best when QC count is moderate to high and spacing is reasonably tight. Intermediate flexibility."
      )
    ))
  }

  if (allow_random_forest) {
    items <- c(items, list(
      htmltools::tags$li(
        htmltools::tags$strong("Random forest (QC-RFSC): "),
        "Best when QC count is high and QC spacing is sufficiently dense. Highest flexibility and highest QC support requirement."
      )
    ))
  }

  caution_block <- htmltools::tags$div(
    style = "margin-top: 8px;",
    warn_span(
      show = total_qcs < 5,
      text = "Flexible methods are not recommended here because the number of QCs is low."
    ),
    warn_span(
      show = allow_random_forest && !spacing_dense_rf,
      text = "Random forest is available by QC count, but current QC spacing does not strongly support it."
    ),
    warn_span(
      show = spacing_sparse,
      text = sprintf(
        "Wide QC spacing detected (max gap = %s). Prefer more stable methods.",
        format(max_gap, digits = 3)
      )
    )
  )

  htmltools::tagList(
    summary_bits,
    rec_block,
    # htmltools::tags$ul(items),
    caution_block
  )
}

# Unavailable correction option description for section 2.1 Choose Correction settings
ui_unavailable_options <- function(df, metab_cols) {
  # If there is only 1 class: class/metabolite impute for samples not available
  # If there is only 1 batch: no batchwise options
  # If there is less than 5 QCs per batch: no batchwise options
  # If there is only 1 batch and less than 5 QCs no RF option
  qc_per_batch <- df |>
    dplyr::group_by(batch) |>
    dplyr::summarise(qc_in_batch = sum(class == "QC"), .groups = "drop")
  num_batches <- length(unique(df$batch))

  sam_df <- df |>
    dplyr::filter(df$class != "QC")
  has_sam_na <- any(is.na(sam_df[, metab_cols]))
  num_classes <- length(unique(sam_df$class))

  unavail_opts <- list()
  if (has_sam_na & (num_classes == 1)) {
    unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Sample Imputation Methods:")
    unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
      icon("circle-xmark", class = "text-danger-emphasis"),
      " class-metabolite median requires more than 1 class."
    )
    unavail_opts[[length(unavail_opts) + 1]] <- tags$br()
    unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
      icon("circle-xmark", class = "text-danger-emphasis"),
      " class-metabolite mean requires more than 1 class."
    )
  }
  if (num_batches == 1) {
    if (any(qc_per_batch$qc_in_batch <= 5)) {
      unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Correction Methods:")
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Random Forest requires more QC samples."
      )
      unavail_opts[[length(unavail_opts) + 1]] <- tags$br()
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise options (Random Forest and LOESS) require more than 1 batch."
      )
    } else {
      unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Correction Methods:")
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise options (Random Forest and LOESS) require more than 1 batch."
      )
    }
  } else {
    if (any(qc_per_batch$qc_in_batch < 5)) {
      unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Correction Methods:")
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise Random Forest requires at least 5 QCs per batch."
      )
      unavail_opts[[length(unavail_opts) + 1]] <- tags$br()
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise LOESS requires at least 5 QCs per batch."
      )
    }
  }

  if (length(unavail_opts) == 0) {
    return(tags$span("All methods available"))
  } else if (length(unavail_opts) == 1) {
    return(unavail_opts[[1]])
  } else {
    return(do.call(tagList, unavail_opts))
  }
}

# ------------------------------------------------------------------------------
# UI RSD summary table
# ------------------------------------------------------------------------------

ui_rsd_stats <- function(compare_to, p, d) {
  rsd_data <- .get_rsd_data_after(compare_to, p, d)
  rsd_results <- .build_rsd_results(rsd_data$df_before, rsd_data$df_after)

  met_stats <- rsd_results$metabolite$stats
  class_stats <- rsd_results$class_metabolite$stats

  met_qc <- met_stats[met_stats$Type == "QC", , drop = FALSE]
  met_samples <- met_stats[met_stats$Type == "Samples", , drop = FALSE]
  class_samples <- class_stats[class_stats$Type == "Samples", , drop = FALSE]

  df <- data.frame(
    Metric = c(
      "Median &Delta; QC RSD",
      "Median &Delta; Metabolite RSD",
      "Median &Delta; Class-Metabolite RSD"
    ),
    Value = c(
      if (nrow(met_qc)) met_qc$med_delta else NA_real_,
      if (nrow(met_samples)) met_samples$med_delta else NA_real_,
      if (nrow(class_samples)) class_samples$med_delta else NA_real_
    ),
    check.names = FALSE
  )

  df$Value <- sprintf("%.2f%%", df$Value)

  change_df <- data.frame(
    s_type = c(
      "QC RSD",
      "Metabolite RSD",
      "Class-Metabolite RSD"
    ),
    increased = c(
      if (nrow(met_qc)) met_qc$pct_increase else NA_real_,
      if (nrow(met_samples)) met_samples$pct_increase else NA_real_,
      if (nrow(class_samples)) class_samples$pct_increase else NA_real_
    ),
    decreased = c(
      if (nrow(met_qc)) met_qc$pct_decrease else NA_real_,
      if (nrow(met_samples)) met_samples$pct_decrease else NA_real_,
      if (nrow(class_samples)) class_samples$pct_decrease else NA_real_
    ),
    check.names = FALSE
  )

  change_df$increased <- sprintf("%.1f%%", change_df$increased)
  change_df$decreased <- sprintf("%.1f%%", change_df$decreased)

  htmltools::tagList(
    htmltools::tags$table(
      style = "border-collapse: collapse; margin-top:10px;",
      htmltools::tags$thead(
        htmltools::tags$tr(
          htmltools::tags$th(
            "Performance Metric",
            style = "padding:4px 12px; text-align:left; border-bottom:1px solid #ccc;"
          ),
          htmltools::tags$th(
            "Value",
            style = "padding:4px 12px; text-align:right; border-bottom:1px solid #ccc;"
          )
        )
      ),
      htmltools::tags$tbody(
        lapply(seq_len(nrow(df)), function(i) {
          htmltools::tags$tr(
            htmltools::tags$td(
              htmltools::HTML(df$Metric[i]),
              style = "padding:4px 12px; text-align:left;"
            ),
            htmltools::tags$td(
              df$Value[i],
              style = "padding:4px 12px; text-align:right;"
            )
          )
        })
      )
    ),
    htmltools::tags$p(),
    htmltools::tags$table(
      style = "border-collapse: collapse; margin-top:10px;",
      htmltools::tags$thead(
        htmltools::tags$tr(
          htmltools::tags$th(
            rsd_data$title,
            style = "padding:4px 12px; text-align:left; border-bottom:1px solid #ccc;"
          ),
          htmltools::tags$th(
            "Increased",
            style = "padding:4px 12px; text-align:right; border-bottom:1px solid #ccc;"
          ),
          htmltools::tags$th(
            "Decreased",
            style = "padding:4px 12px; text-align:right; border-bottom:1px solid #ccc;"
          )
        )
      ),
      htmltools::tags$tbody(
        lapply(seq_len(nrow(change_df)), function(i) {
          htmltools::tags$tr(
            htmltools::tags$td(
              htmltools::HTML(change_df$s_type[i]),
              style = "padding:4px 12px; text-align:left;"
            ),
            htmltools::tags$td(
              change_df$increased[i],
              style = "padding:4px 12px; text-align:right;"
            ),
            htmltools::tags$td(
              change_df$decreased[i],
              style = "padding:4px 12px; text-align:right;"
            )
          )
        })
      )
    )
  )
}
