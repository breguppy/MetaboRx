#' Detect metabolites below a blank/QC threshold
#'
#' Identifies metabolites where the QC mean is less than a user-defined multiple
#' of the blank or processing blank mean. Blank means must be finite and greater
#' than zero to be eligible for threshold comparison.
#'
#' @param df A cleaned non-blank data frame containing metadata columns and
#'   metabolite columns.
#' @param blank_df A data frame containing blank or processing blank rows.
#' @param metab_cols Character vector of metabolite column names to evaluate.
#' @param class_col Name of the sample class column.
#' @param qc_label Label used to identify QC samples.
#' @param threshold Numeric multiplier applied to blank means. A metabolite fails
#'   when `qc_mean < threshold * blank_mean`.
#' @param internal_standard_pattern Regex pattern used to identify internal
#'   standard columns.
#'
#' @return A list containing blank means, QC means, failed metabolites, failed
#'   metabolites excluding internal standards, and a summary table.
#'
#' @keywords internal
#' @noRd
detect_blank_threshold <- function(df,
                                   blank_df,
                                   metab_cols,
                                   class_col = "class",
                                   qc_label = "QC",
                                   threshold = 3,
                                   internal_standard_pattern = "ISTD|ITSD") {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }
  
  if (!is.data.frame(blank_df)) {
    stop("`blank_df` must be a data frame.")
  }
  
  if (!class_col %in% names(df)) {
    stop("`class_col` must exist in `df`.")
  }
  
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || threshold <= 0) {
    stop("`threshold` must be a single positive numeric value.")
  }
  
  metab_cols <- intersect(metab_cols, names(df))
  metab_cols <- intersect(metab_cols, names(blank_df))
  
  if (length(metab_cols) == 0L) {
    return(
      list(
        blank_means = numeric(0),
        qc_means = numeric(0),
        below_blank_threshold = character(0),
        below_blank_threshold_ex_ISTD = character(0),
        threshold_table = data.frame()
      )
    )
  }
  
  if (nrow(blank_df) == 0L) {
    return(
      list(
        blank_means = stats::setNames(rep(NA_real_, length(metab_cols)), metab_cols),
        qc_means = stats::setNames(rep(NA_real_, length(metab_cols)), metab_cols),
        below_blank_threshold = character(0),
        below_blank_threshold_ex_ISTD = character(0),
        threshold_table = data.frame(
          metabolite = metab_cols,
          blank_mean = NA_real_,
          qc_mean = NA_real_,
          threshold_value = NA_real_,
          eligible = FALSE,
          below_blank_threshold = FALSE,
          internal_standard = grepl(
            internal_standard_pattern,
            metab_cols,
            ignore.case = TRUE
          ),
          stringsAsFactors = FALSE
        )
      )
    )
  }
  
  qc_idx <- trimws(as.character(df[[class_col]])) == qc_label
  
  if (!any(qc_idx, na.rm = TRUE)) {
    stop("No QC samples found; cannot compute blank threshold.")
  }
  
  blank_means <- vapply(
    blank_df[, metab_cols, drop = FALSE],
    FUN = function(x) {
      mean(suppressWarnings(as.numeric(x)), na.rm = TRUE)
    },
    FUN.VALUE = numeric(1L)
  )
  
  qc_means <- vapply(
    df[qc_idx, metab_cols, drop = FALSE],
    FUN = function(x) {
      mean(suppressWarnings(as.numeric(x)), na.rm = TRUE)
    },
    FUN.VALUE = numeric(1L)
  )
  
  eligible <- is.finite(blank_means) &
    !is.na(blank_means) &
    blank_means > 0
  
  threshold_value <- threshold * blank_means
  
  failed <- eligible & qc_means < threshold_value
  failed[is.na(failed)] <- FALSE
  
  below_blank_threshold <- names(qc_means)[failed]
  
  internal_standard <- grepl(
    internal_standard_pattern,
    names(qc_means),
    ignore.case = TRUE
  )
  
  below_blank_threshold_ex_ISTD <- below_blank_threshold[
    !grepl(
      internal_standard_pattern,
      below_blank_threshold,
      ignore.case = TRUE
    )
  ]
  
  threshold_table <- data.frame(
    metabolite = names(qc_means),
    blank_mean = unname(blank_means),
    qc_mean = unname(qc_means),
    threshold_value = unname(threshold_value),
    eligible = unname(eligible),
    below_blank_threshold = unname(failed),
    internal_standard = unname(internal_standard),
    stringsAsFactors = FALSE
  )
  
  list(
    blank_means = blank_means,
    qc_means = qc_means,
    below_blank_threshold = below_blank_threshold,
    below_blank_threshold_ex_ISTD = below_blank_threshold_ex_ISTD,
    threshold_table = threshold_table
  )
}

#' Remove metabolites that fail blank threshold detection
#'
#' Removes failed metabolite columns from a data frame. By default, internal
#' standards are protected from removal.
#'
#' @param df A data frame.
#' @param failed_cols Character vector of failed metabolite columns.
#' @param protect_internal_standards Logical. If TRUE, columns matching
#'   `internal_standard_pattern` are not removed.
#' @param internal_standard_pattern Regex pattern used to identify internal
#'   standards.
#'
#' @return A list containing the filtered data frame and removed column names.
#'
#' @keywords internal
#' @noRd
apply_blank_threshold_filter <- function(df,
                                         failed_cols,
                                         protect_internal_standards = TRUE,
                                         internal_standard_pattern = "ISTD|ITSD") {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.")
  }
  
  failed_cols <- intersect(failed_cols, names(df))
  
  if (protect_internal_standards) {
    failed_cols <- failed_cols[
      !grepl(internal_standard_pattern, failed_cols, ignore.case = TRUE)
    ]
  }
  
  filtered_df <- df[, !names(df) %in% failed_cols, drop = FALSE]
  
  list(
    df = filtered_df,
    removed_blank_threshold_cols = failed_cols
  )
}

#' Metabolite/value filtering functions
#'
#' @keywords internal
#' @noRd
filter_by_missing <- function(df, metab_cols, mv_cutoff) {
  if (!"class" %in% names(df)) {
    stop("`df` must contain a 'class' column.")
  }
  
  # get metadata columns
  meta_cols <- setdiff(names(df), metab_cols)
  
  # classes used for group-wise missingness checks
  classes_seen <- sort(unique(df$class))
  
  # compute class-specific missing percentages for each metabolite
  # missing is defined as NA or <= 0
  missing_pct_by_class <- lapply(metab_cols, function(met) {
    vals <- df[[met]]
    
    stats::setNames(
      vapply(classes_seen, function(cl) {
        idx <- which(df$class == cl)
        
        if (length(idx) == 0L) {
          return(100)
        }
        
        x <- vals[idx]
        mean(is.na(x) | x <= 0) * 100
      }, numeric(1L)),
      classes_seen
    )
  })
  names(missing_pct_by_class) <- metab_cols
  
  # remove metabolite if any class exceeds mv_cutoff
  mv_keep_cols <- metab_cols[
    vapply(
      missing_pct_by_class,
      function(x) all(x <= mv_cutoff),
      logical(1L)
    )
  ]
  
  # list columns removed due to missing value %
  mv_removed_cols <- setdiff(metab_cols, mv_keep_cols)
  
  # Get all class-metabolite pairs where all values are missing-like
  # among retained metabolites. Missing-like means NA or <= 0.
  class_metab_all_missing <- if (length(classes_seen) == 0L || length(mv_keep_cols) == 0L) {
    data.frame(
      class = character(0),
      metabolite = character(0),
      n_rows_in_class = integer(0)
    )
  } else {
    out <- vector("list", length(classes_seen))
    names(out) <- classes_seen
    
    for (cl in classes_seen) {
      idx <- which(df$class == cl)
      n_in_class <- length(idx)
      
      if (n_in_class == 0L) {
        miss_all <- rep(TRUE, length(mv_keep_cols))
      } else {
        sub <- df[idx, mv_keep_cols, drop = FALSE]
        miss_all <- vapply(
          sub,
          function(x) all(is.na(x) | x <= 0),
          logical(1L)
        )
      }
      
      mets <- mv_keep_cols[miss_all]
      out[[cl]] <- if (length(mets) == 0L) {
        NULL
      } else {
        data.frame(
          class = rep(cl, length(mets)),
          metabolite = mets,
          n_rows_in_class = rep(n_in_class, length(mets)),
          row.names = NULL
        )
      }
    }
    
    do.call(rbind, Filter(Negate(is.null), out)) %||%
      data.frame(
        class = character(0),
        metabolite = character(0),
        n_rows_in_class = integer(0)
      )
  }
  
  # filter data by metabolite missing value
  df_filtered <- df[, c(meta_cols, mv_keep_cols), drop = FALSE]
  
  # Get retained metabolites that have QC missing-like values
  qc_idx <- which(df_filtered$class == "QC")
  if (length(mv_keep_cols) == 0L || length(qc_idx) == 0L) {
    qc_missing_mets <- character(0)
  } else {
    sub <- df_filtered[qc_idx, mv_keep_cols, drop = FALSE]
    qc_missing_mets <- mv_keep_cols[
      vapply(sub, function(x) any(is.na(x) | x <= 0), logical(1L))
    ]
  }
  
  return(list(
    df = df_filtered,
    mv_cutoff = mv_cutoff,
    mv_removed_cols = mv_removed_cols,
    qc_missing_mets = qc_missing_mets,
    class_metab_all_missing = class_metab_all_missing
  ))
}

remove_imputed_from_corrected <- function(raw_df, corrected_df) {
  # Ensure both data frames are the same shape
  if (!all(dim(raw_df) == dim(corrected_df))) {
    stop("Both data frames must have the same dimensions.")
  }
  
  # Return a new corrected_df with values removed where raw_df is NA
  corrected_df[is.na(raw_df)] <- NA
  return(corrected_df)
}

filter_by_qc_rsd <- function(raw_df,
                             corrected_df,
                             rsd_cutoff,
                             remove_imputed,
                             metadata_cols = c("sample", "batch", "class", "order")) {
  # Compute RSD for corrected_df
  rsd_df <- metabolite_rsd(corrected_df, metadata_cols)
  
  # Identify which metabolites to keep and remove
  keep_metabolites <- rsd_df$Metabolite[!is.na(rsd_df$RSD_QC) &
                                          rsd_df$RSD_QC <= rsd_cutoff]
  
  remove_metabolites <- rsd_df$Metabolite[is.na(rsd_df$RSD_QC) |
                                            rsd_df$RSD_QC > rsd_cutoff]
  
  # Columns to retain in filtered data
  final_cols <- c(metadata_cols, keep_metabolites)
  rsd_filtered_df <- corrected_df[, final_cols, drop = FALSE]

  # Return a list with the filtered data and removed metabolites with and without removing imputed values
  return(
    list(
      df_no_mv = rsd_filtered_df,
      df_mv = rsd_filtered_df,
      rsd_cutoff = rsd_cutoff,
      removed_metabolites_no_mv = remove_metabolites,
      removed_metabolites_mv = remove_metabolites
    )
  )
}

#' Identify metabolites with sample average far from QC average
#'
#' Returns metabolite column names where the average across all non-QC rows
#' differs from the average across QC rows by at least `percent_threshold`
#' percent of the QC average.
#'
#' Percent distance is computed as:
#'   abs(sample_mean - qc_mean) / abs(qc_mean) * 100
#'
#' A pseudocount is added to the denominator to avoid division by zero.
#'
#' @param df data.frame
#'   Input data frame containing metadata columns and metabolite columns.
#' @param metab_cols character or NULL, optional
#'   Metabolite column names. If NULL, uses all columns except
#'   `"sample"`, `"batch"`, `"class"`, and `"order"`.
#' @param class_col character, default "class"
#'   Name of the class column.
#' @param qc_label character, default "QC"
#'   Label used to identify QC rows.
#' @param percent_threshold numeric, default 100
#'   Minimum percent distance from the QC average required for a metabolite to be
#'   returned. A value of 100 means the sample average differs from the QC
#'   average by at least 100 percent of the QC average.
#' @param na_rm logical, default TRUE
#'   Whether to remove missing values when computing averages.
#' @param pseudocount numeric, default 1e-8
#'   Small value added to the absolute QC average denominator to avoid
#'   divide-by-zero.
#' @param return_stats logical, default FALSE
#'   If TRUE, returns a data.frame with the sample average, QC average, and
#'   percent distance for all metabolites. If FALSE, returns only flagged
#'   metabolite names.
#'
#' @return character or data.frame
#'   If `return_stats = FALSE`, returns a character vector of flagged metabolite
#'   names. If `return_stats = TRUE`, returns a data.frame with one row per
#'   metabolite.
#'
#' @examples
#' flagged_metabs <- get_metabs_pct_diff_vs_qc_average(df)
#'
#' stats <- get_metabs_pct_diff_vs_qc_average(
#'   df,
#'   percent_threshold = 75,
#'   return_stats = TRUE
#' )
#'
#' @noRd
get_metabs_pct_diff_vs_qc_average <- function(
    df,
    metab_cols = NULL,
    class_col = "class",
    qc_label = "QC",
    percent_threshold = 100,
    na_rm = TRUE,
    pseudocount = 1e-8,
    return_stats = FALSE
) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  
  if (!class_col %in% names(df)) {
    stop(sprintf("`df` must contain a '%s' column.", class_col))
  }
  
  if (
    !is.numeric(percent_threshold) ||
    length(percent_threshold) != 1L ||
    is.na(percent_threshold) ||
    percent_threshold < 0
  ) {
    stop("`percent_threshold` must be a single numeric value >= 0.")
  }
  
  if (
    !is.numeric(pseudocount) ||
    length(pseudocount) != 1L ||
    is.na(pseudocount) ||
    pseudocount < 0
  ) {
    stop("`pseudocount` must be a single numeric value >= 0.")
  }
  
  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("`na_rm` must be TRUE or FALSE.")
  }
  
  if (!is.logical(return_stats) || length(return_stats) != 1L || is.na(return_stats)) {
    stop("`return_stats` must be TRUE or FALSE.")
  }
  
  if (is.null(metab_cols)) {
    meta_cols <- c("sample", "batch", "class", "order")
    metab_cols <- setdiff(names(df), meta_cols)
  }
  
  if (length(metab_cols) == 0L) {
    stop("No metabolite columns were found.")
  }
  
  missing_metabs <- setdiff(metab_cols, names(df))
  
  if (length(missing_metabs) > 0L) {
    stop(
      sprintf(
        "These `metab_cols` are not in `df`: %s",
        paste(missing_metabs, collapse = ", ")
      )
    )
  }
  
  metab_df <- df[, metab_cols, drop = FALSE]
  
  non_numeric <- metab_cols[!vapply(metab_df, is.numeric, logical(1L))]
  
  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "All metabolite columns must be numeric. Non-numeric columns: %s",
        paste(non_numeric, collapse = ", ")
      )
    )
  }
  
  is_qc <- !is.na(df[[class_col]]) & df[[class_col]] == qc_label
  is_sample <- !is.na(df[[class_col]]) & df[[class_col]] != qc_label
  
  if (!any(is_qc)) {
    stop(sprintf("No rows found where `%s == \"%s\"`.", class_col, qc_label))
  }
  
  if (!any(is_sample)) {
    stop(sprintf("No sample rows found where `%s != \"%s\"`.", class_col, qc_label))
  }
  
  sample_df <- df[is_sample, metab_cols, drop = FALSE]
  qc_df <- df[is_qc, metab_cols, drop = FALSE]
  
  sample_means <- colMeans(sample_df, na.rm = na_rm)
  qc_means <- colMeans(qc_df, na.rm = na_rm)
  
  denominator <- abs(qc_means) + pseudocount
  
  percent_distance <- abs(sample_means - qc_means) / denominator * 100
  
  flagged <- percent_distance >= percent_threshold
  
  if (isTRUE(return_stats)) {
    return(
      data.frame(
        metabolite = metab_cols,
        sample_mean = unname(sample_means),
        qc_mean = unname(qc_means),
        percent_distance_from_qc_average = unname(percent_distance),
        flagged = unname(flagged),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    )
  }
  
  names(percent_distance)[flagged]
}


#' Remove metabolites with sample average far from QC average
#'
#' Removes metabolite columns flagged by
#' `get_metabs_pct_diff_vs_qc_average()`. A metabolite is removed when the
#' average intensity across non-QC samples differs from the average QC intensity
#' by at least `percent_threshold` percent of the QC average.
#'
#' @param df data.frame
#'   Input data frame containing metadata columns and metabolite columns.
#' @param metab_cols character or NULL, optional
#'   Metabolite column names. If NULL, uses all columns except
#'   `"sample"`, `"batch"`, `"class"`, and `"order"`.
#' @param class_col character, default "class"
#'   Name of the class column.
#' @param qc_label character, default "QC"
#'   Label used to identify QC rows.
#' @param percent_threshold numeric, default 100
#'   Minimum percent difference from the QC average required for removal.
#' @param na_rm logical, default TRUE
#'   Whether to remove missing values when computing averages.
#' @param pseudocount numeric, default 1e-8
#'   Small value added to the absolute QC average denominator to avoid
#'   divide-by-zero.
#' @param return_result logical, default FALSE
#'   If TRUE, returns a list containing the filtered data frame, removed
#'   metabolite names, retained metabolite names, and metric table.
#'   If FALSE, returns only the filtered data frame.
#'
#' @return data.frame or list
#'   If `return_result = FALSE`, returns the filtered data frame.
#'   If `return_result = TRUE`, returns a list with:
#'   \describe{
#'     \item{df}{Filtered data frame.}
#'     \item{removed_metabolites}{Metabolite columns removed.}
#'     \item{retained_metabolites}{Metabolite columns retained.}
#'     \item{stats}{Metric table returned by `get_metabs_pct_diff_vs_qc_average()`.}
#'   }
#'
#' @examples
#' filtered_df <- remove_metabs_pct_diff_vs_qc_average(df)
#'
#' result <- remove_metabs_pct_diff_vs_qc_average(
#'   df,
#'   percent_threshold = 100,
#'   return_result = TRUE
#' )
#'
#' @noRd
remove_metabs_pct_diff_vs_qc_average <- function(
    df,
    metab_cols = NULL,
    class_col = "class",
    qc_label = "QC",
    percent_threshold = 100,
    na_rm = TRUE,
    pseudocount = 1e-8,
    return_result = FALSE
) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  
  if (is.null(metab_cols)) {
    meta_cols <- c("sample", "batch", "class", "order")
    metab_cols <- setdiff(names(df), meta_cols)
  }
  
  if (length(metab_cols) == 0L) {
    stop("No metabolite columns were found.")
  }
  
  missing_metabs <- setdiff(metab_cols, names(df))
  
  if (length(missing_metabs) > 0L) {
    stop(
      sprintf(
        "These `metab_cols` are not in `df`: %s",
        paste(missing_metabs, collapse = ", ")
      )
    )
  }
  
  stats <- get_metabs_pct_diff_vs_qc_average(
    df = df,
    metab_cols = metab_cols,
    class_col = class_col,
    qc_label = qc_label,
    percent_threshold = percent_threshold,
    na_rm = na_rm,
    pseudocount = pseudocount,
    return_stats = TRUE
  )
  
  removed_metabolites <- stats$metabolite[stats$flagged]
  retained_metabolites <- setdiff(metab_cols, removed_metabolites)
  
  filtered_df <- df[, setdiff(names(df), removed_metabolites), drop = FALSE]
  
  if (!isTRUE(return_result)) {
    return(filtered_df)
  }
  
  list(
    df = filtered_df,
    removed_metabolites = removed_metabolites,
    retained_metabolites = retained_metabolites,
    stats = stats
  )
}