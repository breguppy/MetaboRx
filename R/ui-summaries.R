#' @keywords internal
#' @noRd

# Metric card for 1.2 Select non-metabolite columns
metric_card <- function(label, value) {
  htmltools::tags$div(
    style = "background:#f8f9fa; padding:10px; border-radius:8px; flex:1;",
    p(style = "font-size:1.4em; font-weight:bold; margin:0;", value),
    h5(style = "margin:0;", label)
  )
}
# Download card used for all download buttons on the first 3 tabs
download_card <- function(title,
                      body,
                      btn) {
  
  header <- shiny::tags$div(
    style = "display:flex; align-items:center; justify-content:space-between; gap:8px;",
    shiny::tags$span(title),
  )
  
  shiny::tags$div(
    class = "card bg-light mb-3",
    style = "margin-top: 10px;",
    shiny::tags$div(class = "card-header", header),
    shiny::tags$div(
      class = "card-body",
      shiny::tags$p(class = "card-text", body),
      btn
    )
  )
}

info_card <- function(title,
                      body,
                      body_tags = NULL,
                      info_title = NULL,
                      info_content = NULL,
                      info_placement = "auto",
                      info_class = "popover-responsive") {
  has_info <- !is.null(info_title) || !is.null(info_content)
  
  header <- shiny::tags$div(
    style = "display:flex; align-items:center; justify-content:space-between; gap:8px;",
    shiny::tags$span(title),
    if (has_info) {
      bslib::popover(
        shiny::tags$button(
          type = "button",
          class = "btn btn-link p-0",
          style = "text-decoration:none;",
          shiny::icon("circle-info")
        ),
        info_content %||% shiny::tags$p(""),
        title = info_title %||% "",
        placement = info_placement,
        options = list(container = "body", customClass = info_class)
      )
    }
  )
  
  shiny::tags$div(
    class = "card border-info mb-3",
    #class = "alert alert-info",
    style = "margin-top: 10px;",
    shiny::tags$div(class = "card-header", 
                    header),
    shiny::tags$div(
      class = "card-body",
      shiny::tags$p(class = "card-text", body),
      body_tags
    )
  )
}
# Reusable warning card generator
warn_card <- function(title,
                      body,
                      body_tags = NULL,
                      info_title = NULL,
                      info_content = NULL,
                      info_placement = "auto",
                      info_class = "popover-responsive") {
  has_info <- !is.null(info_title) || !is.null(info_content)
  
  header <- shiny::tags$div(
    style = "display:flex; align-items:center; justify-content:space-between; gap:8px;",
    shiny::tags$span(title),
    if (has_info) {
      bslib::popover(
        shiny::tags$button(
          type = "button",
          class = "btn btn-link p-0",
          style = "text-decoration:none;",
          shiny::icon("circle-info")
        ),
        info_content %||% shiny::tags$p(""),
        title = info_title %||% "",
        placement = info_placement,
        options = list(
          container = "body",
          customClass = info_class
        )
      )
    }
  )
  
  shiny::tags$div(
    class = "card border-warning mb-3",
    style = "margin-top: 10px;",
    shiny::tags$div(class = "card-header",
                    header),
    shiny::tags$div(
      class = "card-body",
      shiny::tags$p(class = "card-text", body),
      body_tags
    )
  )
}

# Basic info for section 1.2 Select non-metabolite columns
ui_basic_info <- function(cleaned) {
  
  df <- cleaned$df
  replacement_counts <- cleaned$replacement_counts
  non_numeric_cols <- cleaned$non_numeric_cols
  all_missing_zero_qc_cols <- cleaned$all_missing_zero_qc_cols
  duplicate_mets <- cleaned$duplicate_mets
  duplicate_col_names <- cleaned$duplicate_col_names
  blank_df <- cleaned$blank_df
  below_blank_threshold <- cleaned$below_blank_threshold_ex_ISTD
  
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  n_metab    <- length(metab_cols)
  n_missv    <- sum(is.na(df[, metab_cols]))
  n_qcs      <- sum(df$class == "QC")
  n_samp     <- sum(df$class != "QC")
  n_bat      <- dplyr::n_distinct(df$batch)
  n_class    <- dplyr::n_distinct(df$class[df$class != "QC"])
  class_list <- sort(unique(df$class[df$class != "QC"]))
  perc_missv <- round(100 * (n_missv / ((n_samp + n_qcs) * n_metab)), digits = 2)
  
  qc_per_batch <- df %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(qc_in_class = sum(class == "QC"), .groups = "drop")
  
  total_replaced <- sum(
    replacement_counts$non_numeric_replaced +
      replacement_counts$zero_replaced
  )
  
  class_badges <- tags$div(
    style = "display: flex; flex-wrap: wrap; gap: 8px; margin-top: 5px;",
    lapply(class_list, function(cls) {
      tags$span(
        style = "background-color: #e9ecef; padding: 5px 10px; border-radius: 12px;",
        as.character(cls)
      )
    })
  )
  
  # ---------- Warning box 1: replaced values ----------
  replaced_card <- NULL
  if (total_replaced > 0) {
    replaced_card <- info_card(
      title = "Replaced non-numeric or zeros values in metabolite columns",
      body  = paste0(
        total_replaced,
        " values were converted to missing (NA) prior to processing."
      ),
    )
  }
  
  # ---------- Warning box 2: removed metabolite columns ----------
  nonnum_card <- NULL
  
  removed_non_numeric <- sort(unique(non_numeric_cols))
  removed_all_zero_qc <- sort(unique(all_missing_zero_qc_cols))
  
  if (
    length(removed_non_numeric) > 0 ||
    length(removed_all_zero_qc) > 0
  ) {
    
    section_tag <- function(title, values) {
      if (length(values) == 0) return(NULL)
      
      tags$div(
        tags$p(
          style = "font-weight: 600; margin-top: 8px; margin-bottom: 4px;",
          title
        ),
        tags$p(
          style = "margin-bottom: 6px;",
          paste(values, collapse = ", ")
        )
      )
    }
    
    nonnum_card <- warn_card(
      title = "Removed metabolite columns",
      body  = "The following metabolite columns were removed prior to processing:",
      body_tags = tags$div(
        section_tag(
          "Non-numerical columns:",
          removed_non_numeric
        ),
        section_tag(
          "All values missing or zero for QC samples:",
          removed_all_zero_qc
        )
      )
    )
  }
  
  # ---------- Warning box 3: duplicate column names ----------
  duplicate_columns <- NULL
  if (!is.null(duplicate_col_names) && length(duplicate_col_names) > 0) {
    
    duplicate_columns <- warn_card(
      title = "Duplicate column names",
      body  = "The follow column names appear more than once in your dataset. We have appended '_1', '_2', etc. to subsequent duplicates.",
      body_tags = tags$ul(
        style = "margin-bottom: 0;",
        lapply(sort(duplicate_col_names), tags$li)
      )
    )
  }
  
  # ---------- Warning box 4: duplicate metabolites ----------
  duplicate_card <- NULL
  if (!is.null(duplicate_mets) && nrow(duplicate_mets) > 0) {
    
    dup_badges <- tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px;",
      lapply(seq_len(nrow(duplicate_mets)), function(i) {
        pair <- duplicate_mets[i, , drop = FALSE]
        tags$span(
          style = paste(
            "background-color: #ff8989;",
            "border: 1px solid #ff8989;",
            "padding: 4px 8px;",
            "border-radius: 12px;",
            "font-size: 0.85rem;"
          ),
          sprintf("%s \u2248 %s", pair$col1, pair$col2)
        )
      })
    )
    
    duplicate_card <- warn_card(
      title = "Potential duplicate metabolites",
      body  = sprintf(
        "%d column pairs appear equal or nearly equal based on non-missing values.",
        nrow(duplicate_mets)
      ),
      body_tags = dup_badges,
    )
  }
  # ---------------------------------------------------
  # MAIN UI SECTION
  # ---------------------------------------------------
  tagList(
    replaced_card,
    nonnum_card,
    duplicate_columns,
    duplicate_card,
    
    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 20px; margin-top: 10px;",
      
      # left section
      tags$div(
        style = "display: grid; grid-template-columns: repeat(1, 1fr); gap: 20px;",
        
        # metrics grid
        tags$div(
          style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin-top: 15px;",
          metric_card("Metabolite Columns", n_metab),
          metric_card("Missing Values", paste0(n_missv, " (", perc_missv, "%)")),
          metric_card("QC Samples", n_qcs),
          metric_card("Samples", n_samp),
          metric_card("Batches", n_bat),
          metric_card("Classes", n_class)
        ),
        
        # class list badges
        tags$div(
          style = "flex: 1; min-width: 250px;",
          tags$h5("Unique Classes"),
          class_badges
        )
      ),
      
      # right section
      tags$div(
        style = "flex: 1; min-width: 250px;",
        tags$h5("Number of QC Samples per Batch"),
        tags$table(
          class = "table table-bordered table-sm",
          tags$thead(
            tags$tr(tags$th("Batch"), tags$th("QCs in Batch"))
          ),
          tags$tbody(
            lapply(seq_len(nrow(qc_per_batch)), function(i) {
              tags$tr(
                tags$td(as.character(qc_per_batch$batch[i])),
                tags$td(qc_per_batch$qc_in_class[i])
              )
            })
          )
        )
      )
    )
  )
}

#' Blank threshold filtering info card
#'
#' @param blank_threshold_result Result from detect_blank_threshold().
#' @param blank_df Data frame containing blank or processing blank samples.
#' @param threshold Numeric blank threshold multiplier.
#' @param remove_blank_threshold_cols Logical indicating whether failed columns
#'   are removed.
#' @param removed_blank_threshold_cols Character vector of removed metabolite
#'   columns.
#'
#' @return A Shiny tag object or NULL.
#'
#' @keywords internal
#' @noRd
ui_blank_threshold_info <- function(blank_threshold_result,
                                    blank_df,
                                    threshold,
                                    remove_blank_threshold_cols = FALSE,
                                    removed_blank_threshold_cols = character(0)) {
  n_blanks <- if (!is.null(blank_df) && nrow(blank_df) > 0L) {
    nrow(blank_df)
  } else {
    0L
  }
  
  if (n_blanks == 0L || is.null(blank_threshold_result)) {
    return(NULL)
  }
  
  below_blank_threshold <- unique(stats::na.omit(
    as.character(blank_threshold_result$below_blank_threshold_ex_ISTD)
  ))
  
  blank_body <- sprintf(
    "%d blank/processing blank sample(s) detected and excluded from downstream processing.",
    n_blanks
  )
  
  threshold_status <- if (length(below_blank_threshold) > 0L) {
    shiny::tags$div(
      shiny::tags$p(
        style = "font-weight: 600; margin-top: 8px; margin-bottom: 6px;",
        sprintf(
          "%d metabolite(s) failed the %.1fx blank-average threshold for QC samples, excluding internal standards:",
          length(below_blank_threshold),
          threshold
        )
      ),
      shiny::tags$ul(
        style = "margin-bottom: 0;",
        lapply(sort(below_blank_threshold), shiny::tags$li)
      )
    )
  } else {
    shiny::tags$p(
      style = "font-weight: 600; margin-top: 8px; margin-bottom: 0;",
      sprintf(
        "All metabolites have QC average above %.1f times the average of blanks.",
        threshold
      )
    )
  }
  
  removal_status <- if (isTRUE(remove_blank_threshold_cols)) {
    if (length(removed_blank_threshold_cols) > 0L) {
      shiny::tags$p(
        style = "margin-top: 8px; margin-bottom: 0;",
        sprintf(
          "%d metabolite column(s) were removed before missing-value filtering.",
          length(removed_blank_threshold_cols)
        )
      )
    } else {
      shiny::tags$p(
        style = "margin-top: 8px; margin-bottom: 0;",
        "Blank-threshold removal is enabled, but no metabolite columns were removed."
      )
    }
  } else {
    shiny::tags$p(
      style = "margin-top: 8px; margin-bottom: 0;",
      "Blank-threshold removal is disabled. Failed metabolites are flagged but retained."
    )
  }
  
  warn_card(
    title = "Blank threshold filtering",
    body = blank_body,
    body_tags = shiny::tags$div(
      threshold_status,
      removal_status
    )
  )
}

# Filter info for section 1.4 Filter Raw Data
ui_filter_info <- function(fd, mv_cutoff) {
  mv_removed <- fd$mv_removed_cols
  qc_missing_mets <- fd$qc_missing_mets
  class_metab_all_missing <- fd$class_metab_all_missing
  df <- fd$df
  
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  n_metab    <- length(metab_cols)
  n_missv    <- sum(is.na(df[, metab_cols]))
  n_qcs      <- sum(df$class == "QC")
  n_samp     <- sum(df$class != "QC")
  perc_missv <- round(100 * (n_missv / ((n_samp + n_qcs) * n_metab)), digits = 2)
  
  
  left_col <- if (length(mv_removed) == 0) {
    tags$div(style = "flex: 1; padding-right: 10px;",
             tags$span(
               style = "color:darkgreen;font-weight:bold;",
               paste0(
                 "No metabolites with missing value percentage above ",
                 mv_cutoff,
                 "%."
               )
             ))
  } else {
    tags$div(
      style = "flex: 1; padding-right: 10px;",
      tags$span(
        style = "color:darkorange;font-weight:bold;",
        paste0(
          length(mv_removed),
          " metabolite(s) removed based on missing value percentage above ",
          mv_cutoff,
          "%."
        )
      ),
      tags$ul(lapply(mv_removed, tags$li))
    )
  }
  
  right_col <- if(length(qc_missing_mets) == 0) {
    tags$div(
      class = "alert alert-success",
      style = "margin-bottom: 10px;",
      #tags$span(style = "color:darkgreen; font-weight:bold;",
       tags$strong("No metabolites have missing values in QC samples after filtering."))
  } else {
    tags$div(
      class = "alert alert-warning",
      style = "margin-bottom: 10px;",
      #tags$span(style = "color:darkorange; font-weight:bold;",
      tags$strong(          paste0(length(qc_missing_mets),
                       " metabolite(s) with at least one QC missing value after filtering.")),
      tags$ul(lapply(qc_missing_mets, tags$li)))
  }
  right_col1 <- tags$div(style = "flex: 1; min-width: 250px;",
                        right_col,
                        tags$div(
                          style = "display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px; margin-top: 15px;",
                          metric_card("Metabolite Columns", n_metab),
                          metric_card("Missing Values", paste0(n_missv, " (", perc_missv, "%)")))
  )
  summary_row <- tags$div(
    style = "display:flex; gap:16px; align-items:flex-start;",
    left_col, right_col1
  )
  
  has_all_missing <- !is.null(class_metab_all_missing) &&
    is.data.frame(class_metab_all_missing) &&
    nrow(class_metab_all_missing) > 0L
  
  all_missing_card <- NULL
  if (has_all_missing) {
    # Create bullet list like: "QC — MetaboliteA"
    pair_items <- apply(
      class_metab_all_missing[, c("class", "metabolite"), drop = FALSE],
      1,
      function(r) paste0(r[[1]], " — ", r[[2]])
    )
    
    all_missing_card <- warn_card(
      title = "All-missing class/metabolite combinations detected",
      body  = paste0(
        "The following class–metabolite pairs have all values missing. ",
        "These values will remain missing if you choose a class-metabolite imputation method."
      ),
      body_tags = shiny::tags$ul(lapply(pair_items, shiny::tags$li))
    )
  }
  
  shiny::tagList(
    summary_row,
    all_missing_card
  )
}


.make_correlation_card <- function(high_corr_mets, range, card_title) {
  cor_badges <- NULL
  card_body <- tags$div(
    style = "flex: 1; padding-right: 10px;",
    tags$span(style = "color:darkgreen;font-weight:bold;", "No metabolite pairs within this correlation range.")
  )
  if (!is.null(high_corr_mets) && nrow(high_corr_mets) > 0) {
    cor_badges <- tags$div(style = "display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px;", lapply(seq_len(nrow(high_corr_mets)), function(i) {
      pair <- high_corr_mets[i, , drop = FALSE]
      tags$span(
        style = paste(
          "background-color: #fff3cd;",
          "border: 1px solid #ffeeba;",
          "padding: 4px 8px;",
          "border-radius: 12px;",
          "font-size: 0.85rem;"
        ),
        sprintf("%s \u221D %s", pair$col1, pair$col2)
      )
    }))
    card_body <- sprintf(
      "%d column pairs have correlation (Pearson's r) within the selected range: %.3f - %.3f.",
      nrow(high_corr_mets),
      range[1],
      range[2]
    )
  }
  correlated_card <- info_card(title = card_title,
                               body  = card_body,
                               body_tags = cor_badges)
}
ui_corr_range_info <- function(all_corr, range) {
  raw_high_corr_mets       <- filter_correlation_pairs_by_range(all_corr$raw, range)
  corrected_high_corr_mets <- filter_correlation_pairs_by_range(all_corr$corrected, range)
  
  raw_corr_card       <- .make_correlation_card(raw_high_corr_mets, range, "Raw Metabolite Correlations")
  corrected_corr_card <- .make_correlation_card(corrected_high_corr_mets, range, "Corrected Metabolite Correlations")
  
  transformed_block <- if (isTRUE(all_corr$transformed_included)) {
    transformed_high_corr_mets <- filter_correlation_pairs_by_range(all_corr$transformed, range)
    .make_correlation_card(
      transformed_high_corr_mets,
      range,
      "Transformed and Corrected Metabolite Correlations"
    )
  } else {
    NULL
  }
  
  htmltools::tagList(
    raw_corr_card,
    corrected_corr_card,
    transformed_block
  )
}

# Post-correction filtering info for section 2.2 Post-Correction Filtering
ui_postcor_filter_info <- function(filtered_corrected_result,
                                   remove_imputed,
                                   rsd_cutoff,
                                   post_cor_filter,
                                   remove_qc_average_pct_filter) {
  if (isTRUE(remove_imputed)) {
    removed <- filtered_corrected_result$removed_metabolites_mv
    df <- filtered_corrected_result$df_mv
    metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
    n_missv    <- sum(is.na(df[, metab_cols]))
    imputed_removed_ui <- htmltools::tagList(metric_card(
      "Imputed values are removed after correction",
      n_missv
    ),
    tags$br())
  } else {
    removed <- filtered_corrected_result$removed_metabolites_no_mv
    df <- filtered_corrected_result$df_no_mv
    imputed_removed_ui <- NULL
  }
  
  flagged <- filtered_corrected_result$flagged_mets
  
  n_removed <- length(removed)
  
  # get ISTD/ITSD metabolites
  is_istd <- grepl("ISTD|ITSD", removed, ignore.case = TRUE)
  istd_names <- removed[is_istd]
  n_istd <- length(istd_names)
  
  met_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  total <- n_removed + length(met_cols)
  
  pct_below <- if (total > 0) {
    round((length(met_cols) / total) * 100, digits = 1)
  } else {
    NA_real_
  }
  
  # optional warning banner for internal standards failing RSD filter
  warning_ui <- NULL
  if (n_istd > 0) {
    warning_ui <- tags$div(
      class = "alert alert-danger",
      style = "margin-bottom: 10px;",
      tags$strong(
        paste0(
          n_istd,
          " internal standard(s) with QC RSD above ",
          rsd_cutoff,
          "%:"
        )
      ),
      tags$ul(
        lapply(istd_names, tags$li)
      )
    )
  }
  
  removal_status <- if (isTRUE(remove_qc_average_pct_filter)) {
    if (length(flagged) > 0L) {
      shiny::tags$p(
        style = "margin-top: 8px; margin-bottom: 0;",
        sprintf(
          "%d metabolite column(s) were removed before RSD filtering.",
          length(flagged)
        )
      )
    } else {
      shiny::tags$p(
        style = "margin-top: 8px; margin-bottom: 0;",
        "Differs from QC average removal is enabled, but no metabolite columns were removed."
      )
    }
  } else {
    shiny::tags$p(
      style = "margin-top: 8px; margin-bottom: 0;",
      "Differs from QC average removal is disabled. Failed metabolites are flagged but retained."
    )
  }
  # optional warning banner for metabolites not within 2-fold of QC
  flagged_ui <- NULL
  if (!is.null(flagged) && length(flagged) > 0) {
    flagged_ui <- tags$div(
      class = "alert alert-warning",
      style = "margin-bottom: 10px;",
      tags$strong(
        paste0(
          length(flagged),
          " metabolite(s) where the sample average differs from the QC average by at least ",
          filtered_corrected_result$percent_threshold,
          "%."
        )
      ),
      tags$ul(
        lapply(flagged, tags$li)
      ),
      removal_status
    )
  }
  
  if (isFALSE(post_cor_filter)) {
    ui <- list(
      imputed_removed_ui,
      warning_ui,
      flagged_ui,
      metric_card(
        paste0("Metabolites with QC RSD at or below ", rsd_cutoff, "%"),
        paste0(pct_below, "%")
      ),
      tags$span(
        style = "color: darkorange; font-weight: bold;",
        paste0(
          n_removed,
          " metabolite(s) removed based on QC RSD above ",
          rsd_cutoff,
          "%"
        )
      ),
      tags$br(),
      tags$ul(
        lapply(removed, tags$li)
      )
    )
  } else {
    ui <- list(
      imputed_removed_ui,
      warning_ui,
      flagged_ui,
      tags$span(
        style = "color: darkgreen; font-weight: bold;",
        "Metabolites are not filtered by QC RSD."
      )
    )
  }
  
  do.call(tagList, ui)
}

#' Candidate extreme values summary table UI
#'
#' @param detect_result Result returned by `detect_hotelling_nonqc_dual_z()`.
#' @param top_n Number of candidate extreme values to display.
#' @param sample_col Name of the sample column in the extreme values table.
#' @param class_col Name of the class column in the extreme values table.
#' @param digits_z Number of digits used for z-score display.
#' @param digits_T2 Number of digits used for Mahalanobis distance display.
#'
#' @return A Shiny tag list containing metric cards and the outlier table.
#'
#' @keywords internal
#' @noRd
ui_outliers_table <- function(detect_result,
                              top_n     = 10L,
                              sample_col = "sample",
                              class_col  = "class",
                              digits_z   = 2L,
                              digits_T2  = 2L) {
  if (is.null(detect_result)) {
    stop("detect_result is NULL.")
  }
  
  ev <- detect_result$extreme_values
  dres <- detect_result$data
  
  if (is.null(ev)) {
    stop("detect_result$extreme_values is NULL. Did you pass the correct object?")
  }
  
  if (is.null(dres)) {
    stop("detect_result$data is NULL. Did you pass the correct object?")
  }
  
  n_outlier_samples <- sum(dres$is_outlier_sample, na.rm = TRUE)
  n_extreme_values <- nrow(ev)
  
  cards <- shiny::div(
    style = paste(
      "display:flex;",
      "gap:10px;",
      "margin-bottom:10px;",
      "flex-wrap:wrap;"
    ),
    metric_card("Samples outside the Mahalanobis 95% limit", n_outlier_samples),
    metric_card("Potential extreme metabolite values", n_extreme_values)
  )
  
  if (nrow(ev) == 0L) {
    return(
      shiny::tagList(
        cards,
        shiny::tags$em("No extreme metabolite values detected.")
      )
    )
  }
  
  required_cols <- c(
    sample_col,
    class_col,
    "metabolite",
    "z_global",
    "abs_z_global",
    "z_class",
    "abs_z_class",
    "T2"
  )
  
  missing_cols <- setdiff(required_cols, names(ev))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in extreme_values: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  top_n <- as.integer(top_n)
  
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L) {
    top_n <- 10L
  }
  
  ev_sorted <- ev[order(-ev$abs_z_global, -ev$abs_z_class, -ev$T2), , drop = FALSE]
  ev_top <- head(ev_sorted, top_n)
  
  z_g_fmt <- formatC(ev_top$z_global, format = "f", digits = digits_z)
  z_c_fmt <- formatC(ev_top$z_class, format = "f", digits = digits_z)
  T2_fmt <- formatC(ev_top$T2, format = "f", digits = digits_T2)
  
  rows <- lapply(seq_len(nrow(ev_top)), function(i) {
    shiny::tags$tr(
      shiny::tags$td(ev_top[[sample_col]][i]),
      shiny::tags$td(ev_top[[class_col]][i]),
      shiny::tags$td(ev_top$metabolite[i]),
      shiny::tags$td(z_g_fmt[i]),
      shiny::tags$td(z_c_fmt[i]),
      shiny::tags$td(T2_fmt[i])
    )
  })
  
  table_tag <- shiny::tags$table(
    class = "table table-striped table-condensed table-hover",
    shiny::tags$thead(
      shiny::tags$tr(
        shiny::tags$th("Sample"),
        shiny::tags$th("Class"),
        shiny::tags$th("Metabolite"),
        shiny::tags$th("Global z-score"),
        shiny::tags$th("Class z-score"),
        shiny::tags$th("Mahalanobis^2")
      )
    ),
    shiny::tags$tbody(rows)
  )
  
  shiny::tagList(
    cards,
    shiny::tags$p(
      sprintf("Top %s potential extreme values are listed below. ", top_n),
      "The full list of potential extreme values ",
      sprintf("'extreme_values_%s.xlsx' ", Sys.Date()),
      "is available for download."
    ),
    shiny::tags$div(
      style = "overflow-x:auto; width:100%;",
      table_tag
    )
  )
}
