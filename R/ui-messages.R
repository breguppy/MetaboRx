#' @keywords internal
#' @noRd

# Column-warning text for section 1.2 Select non-metabolite columns
ui_column_warning <- function(data, selected) {
  warnings <- list()
  
  if (any(selected == "")) {
    warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkorange; font-weight:bold;",
                                                  icon("exclamation-triangle"),
                                                  " Please select all four columns.")
  } else if (length(unique(selected)) < 4) {
    warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkred; font-weight:bold;",
                                                  icon("exclamation-triangle"),
                                                  " Each selected column must be unique.")
  }
  
  if (!any(selected == "") && length(unique(selected)) == 4) {
    samp_vec  <- data[[selected[1]]]
    order_vec <- data[[selected[4]]]
    
    if (anyDuplicated(samp_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
                                                    icon("exclamation-triangle"),
                                                    " Duplicate sample names detected!")
    }
    
    if (anyDuplicated(order_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
                                                    icon("exclamation-triangle"),
                                                    " Duplicate order values detected!")
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
  n_missv = sum(is.na(df[qc_idx, metab_cols]))
  
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
  
  total_qcs <- sum(df[[class_col]] == qc_label, na.rm = TRUE)
  
  # Compute QC spacing (assumes order_col exists and is usable at this point)
  qc_gap_stats <- NULL
  if (order_col %in% names(df)) {
    ord <- df[[order_col]]
    is_qc <- df[[class_col]] == qc_label
    
    keep <- !is.na(ord) & !is.na(is_qc)
    ord <- ord[keep]
    is_qc <- is_qc[keep]
    
    if (!is.numeric(ord)) {
      ord_num <- suppressWarnings(as.numeric(ord))
      if (!all(is.na(ord_num))) ord <- ord_num
    }
    
    if (is.numeric(ord) && sum(is_qc) >= 2L) {
      qc_orders <- sort(ord[is_qc])
      gaps <- diff(qc_orders)
      qc_gap_stats <- list(
        max_gap = max(gaps, na.rm = TRUE),
        median_gap = stats::median(gaps, na.rm = TRUE)
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
    "QC spacing (injection order): unavailable (need ≥2 QCs with valid order)"
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
  
  # Helper to render conditional warning text
  warn_span <- function(show, text) {
    if (!isTRUE(show)) return(NULL)
    htmltools::tags$span(
      icon("exclamation-triangle", class = "text-danger-emphasis"),
      htmltools::tags$span(style = "margin-left: 6px;", text),
      style = "display:block; margin-top:4px;"
    )
  }
  
  # Decide overfit-risk flags (data-driven when possible)
  max_gap <- qc_gap_stats$max_gap %||% NA_real_
  
  # LOESS polynomial warning triggers
  loess2_warn_low_qc <- total_qcs < 9
  loess2_warn_gap <- is.finite(max_gap) && max_gap > 15
  
  # RF warning triggers (more conservative)
  rf_warn_low_qc <- total_qcs < 12
  rf_warn_gap <- is.finite(max_gap) && max_gap > 10
  
  # Base choices always available
  items <- list(
    htmltools::tags$li(
      htmltools::tags$strong("Local constant: "),
      "Use when the QC drift trend is flat, dominated by noise, or shows no consistent pattern."
    ),
    htmltools::tags$li(
      htmltools::tags$strong("Local linear: "),
      "Use when the QC drift trend is a gradual increase or decrease (approximately monotone)."
    )
  )
  
  # Add local polynomial when allowed
  if (total_qcs >= 5 && total_qcs <= 8) {
    items <- c(items, list(
      htmltools::tags$li(
        htmltools::tags$strong("Local polynomial (QC-RLSC): "),
        "Use when the QC drift trend is a smooth curve (nonlinear but smooth).",
        warn_span(
          show = loess2_warn_low_qc || loess2_warn_gap,
          text = paste(
            "Higher overfit risk with sparse QCs.",
            if (loess2_warn_low_qc) "With <9 QCs, polynomial fits can be unstable." else "",
            if (loess2_warn_gap) sprintf("Your max QC gap is %s (>15).", format(max_gap, digits = 3)) else "",
            "If you see oscillations or worse non-QC variability, prefer degree 1."
          )
        )
      )
    ))
  }
  
  # Add polynomial + RF when allowed
  if (total_qcs > 8 && total_qcs <= 15) {
    items <- c(items, list(
      htmltools::tags$li(
        htmltools::tags$strong("Local polynomial (QC-RLSC): "),
        "Use when drift is smooth and curved. Prefer degree 1 if drift is mostly linear.",
        warn_span(
          show = loess2_warn_gap,
          text = paste(
            if (loess2_warn_gap) sprintf("Your max QC gap is %s (>15), which can make polynomial LOESS unstable.", format(max_gap, digits = 3)) else "",
            "If the correction curve looks wiggly, increase smoothing or use degree 1."
          )
        )
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Random forest (QC-RFSC): "),
        "Use when drift is irregular, has local disruptions, or shows abrupt changes that a smooth curve cannot capture.",
        warn_span(
          show = rf_warn_low_qc || rf_warn_gap,
          text = paste(
            "Random forest is high-flexibility and can overfit with limited QC support.",
            if (rf_warn_low_qc) "With <12 QCs, overfit risk is elevated." else "",
            if (rf_warn_gap) sprintf("Your max QC gap is %s (>10).", format(max_gap, digits = 3)) else "",
            "Prefer LOESS (degree 1/2) unless LOESS fails to reduce QC drift/RSD."
          )
        )
      )
    ))
  }
  
  if (total_qcs > 15) {
    items <- c(items, list(
      htmltools::tags$li(
        htmltools::tags$strong("Local polynomial (QC-RLSC): "),
        "Use when drift is smooth and curved.",
        warn_span(
          show = loess2_warn_gap,
          text = if (loess2_warn_gap) {
            sprintf("Your max QC gap is %s (>15). Large gaps can make LOESS less reliable between QC anchors.", format(max_gap, digits = 3))
          } else {
            NULL
          }
        )
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Random forest (QC-RFSC): "),
        "Use when drift is irregular or has abrupt changes. Often a good option when QC frequency is high.",
        warn_span(
          show = rf_warn_gap,
          text = if (rf_warn_gap) {
            sprintf("Your max QC gap is %s (>10). Wide QC spacing increases RF overfit risk.", format(max_gap, digits = 3))
          } else {
            NULL
          }
        )
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Rule of thumb: "),
        "If both LOESS polynomial and random forest are available, try both and compare QC drift/RSD reduction and non-QC stability."
      )
    ))
  }
  
  htmltools::tagList(
    summary_bits,
    htmltools::tags$ul(items)
  )
}

# Unavailable correction option description for section 2.1 Choose Correction settings
ui_unavailable_options <- function(df, metab_cols) {
  # If there is only 1 class: class/metabolite impute for samples not available
  # If there is only 1 batch: no batchwise options
  # If there is less than 5 QCs per batch: no batchwise options
  # If there is only 1 batch and less than 5 QCs no RF option
  qc_per_batch <- df %>%
    group_by(batch) %>%
    summarise(qc_in_batch = sum(class == "QC"), .groups = "drop")
  num_batches <- length(unique(df$batch))
  
  sam_df <- df %>% filter(df$class != "QC")
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

ui_rsd_stats <- function(compare_to, p, d) {
  df_before <- d$filtered$df
  # Determine df_after based on rsd_compare selected by user.
  if (compare_to == "filtered_cor_data") {
    if (isTRUE(p$remove_imputed)) {
      df_after <- d$filtered_corrected$df_mv
    } else {
      df_after <- d$filtered_corrected$df_no_mv
    }
  } else {
    if (isTRUE(p$remove_imputed)) {
      df_after <- d$transformed$df_mv
    } else {
      df_after <- d$transformed$df_no_mv
    }
  }
  
  met_rsdBefore <- metabolite_rsd(df_before)
  met_rsdAfter <- metabolite_rsd(df_after)
  
  class_rsdBefore <- class_metabolite_rsd(df_before)
  class_rsdAfter <- class_metabolite_rsd(df_after)
  
  rsd_met_stats <- delta_rsd_stats(met_rsdBefore, met_rsdAfter)
  rsd_class_stats <- delta_rsd_stats(class_rsdBefore, class_rsdAfter)
  
  df <- data.frame(
    Metric = c("Median &Delta; QC RSD", 
               "Median &Delta; Metabolite RSD",
               "Median &Delta; Class-Metabolite RSD"),
    Value  = c(rsd_met_stats$med_delta_qc,
               rsd_met_stats$med_delta_sample,
               rsd_class_stats$med_delta_sample)
  )
  
  df$Value <- sprintf("%.2f%%", df$Value)
  
  change_df <- data.frame(
    s_type = c("QC RSD", 
               "Metabolite RSD",
               "Class-Metabolite RSD"),
    increased  = c(rsd_met_stats$pct_increase_qc,
                   rsd_met_stats$pct_increase_sample,
                   rsd_class_stats$pct_increase_sample),
    decreased  = c(rsd_met_stats$pct_decrease_qc,
                   rsd_met_stats$pct_decrease_sample,
                   rsd_class_stats$pct_decrease_sample)
  )
  
  change_df$increased <- sprintf("%.1f%%", change_df$increased)
  change_df$decreased <- sprintf("%.1f%%", change_df$decreased)
  
  htmltools::tagList(
    
    htmltools::tags$table(
      style = "border-collapse: collapse; margin-top:10px;",
      htmltools::tags$thead(
        htmltools::tags$tr(
          htmltools::tags$th("Performance Metric",  style="padding:4px 12px; text-align:left; border-bottom:1px solid #ccc;"),
          htmltools::tags$th("Value",   style="padding:4px 12px; text-align:right; border-bottom:1px solid #ccc;")
        )
      ),
      htmltools::tags$tbody(
        lapply(seq_len(nrow(df)), function(i) {
          htmltools::tags$tr(
            htmltools::tags$td(HTML(df$Metric[i]), style="padding:4px 12px; text-align:left;"),
            htmltools::tags$td(df$Value[i],  style="padding:4px 12px; text-align:right;")
          )
        })
      )
    ),
    htmltools::tags$p(),
    htmltools::tags$table(
      style = "border-collapse: collapse; margin-top:10px;",
      htmltools::tags$thead(
        htmltools::tags$tr(
          htmltools::tags$th("Post-correction Changes",  style="padding:4px 12px; text-align:left; border-bottom:1px solid #ccc;"),
          htmltools::tags$th("Increased",   style="padding:4px 12px; text-align:right; border-bottom:1px solid #ccc;"),
          htmltools::tags$th("Decreased",   style="padding:4px 12px; text-align:right; border-bottom:1px solid #ccc;")
        )
      ),
      htmltools::tags$tbody(
        lapply(seq_len(nrow(change_df)), function(i) {
          htmltools::tags$tr(
            htmltools::tags$td(HTML(change_df$s_type[i]), style="padding:4px 12px; text-align:left;"),
            htmltools::tags$td(change_df$increased[i],  style="padding:4px 12px; text-align:right;"),
            htmltools::tags$td(change_df$decreased[i],  style="padding:4px 12px; text-align:right;")
          )
        })
      )
    )
  )
}
