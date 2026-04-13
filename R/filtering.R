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
  # corrected_df will have no missing values
  # mv_corrected_df will have imputed values removed.
  mv_corrected_df <- corrected_df
  if (remove_imputed) {
    # Ensure both data frames are the same shape
    if (!all(dim(raw_df) == dim(mv_corrected_df))) {
      stop("Both data frames must have the same dimensions.")
    }
    
    # Return a new corrected_df with values removed where raw_df is NA
    mv_corrected_df[is.na(raw_df)] <- NA
  }

  
  # Compute RSD for corrected_df and mv_corrected_df
  rsd_no_mv_df <- metabolite_rsd(corrected_df, metadata_cols)
  rsd_mv_df <- metabolite_rsd(mv_corrected_df, metadata_cols)
  
  # Identify which metabolites to keep and remove
  keep_metabolites_no_mv <- rsd_no_mv_df$Metabolite[!is.na(rsd_no_mv_df$RSD_QC) &
                                          rsd_no_mv_df$RSD_QC <= rsd_cutoff]
  keep_metabolites_mv <- rsd_mv_df$Metabolite[!is.na(rsd_mv_df$RSD_QC) &
                                                      rsd_mv_df$RSD_QC <= rsd_cutoff]
  
  remove_metabolites_no_mv <- rsd_no_mv_df$Metabolite[is.na(rsd_no_mv_df$RSD_QC) |
                                            rsd_no_mv_df$RSD_QC > rsd_cutoff]
  remove_metabolites_mv <- rsd_mv_df$Metabolite[is.na(rsd_mv_df$RSD_QC) |
                                                        rsd_mv_df$RSD_QC > rsd_cutoff]
  
  # Columns to retain in filtered data
  final_cols_no_mv <- c(metadata_cols, keep_metabolites_no_mv)
  final_cols_mv <- c(metadata_cols, keep_metabolites_mv)
  rsd_filtered_df <- corrected_df[, final_cols_no_mv, drop = FALSE]
  rsd_mv_filtered_df <- mv_corrected_df[, final_cols_mv, drop = FALSE]

  # Return a list with the filtered data and removed metabolites with and without removing imputed values
  return(
    list(
      df_no_mv = rsd_filtered_df,
      df_mv = rsd_mv_filtered_df,
      rsd_cutoff = rsd_cutoff,
      removed_metabolites_no_mv = remove_metabolites_no_mv,
      removed_metabolites_mv = remove_metabolites_mv
    )
  )
}

#' Identify metabolites with >= 2-fold difference between sample and QC averages
#'
#' Returns the metabolite column names where the average across all non-QC rows
#' is at least `fold_threshold` times different from the average across QC rows.
#'
#' Fold-difference is computed symmetrically as:
#'   max(sample_mean / qc_mean, qc_mean / sample_mean)
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
#' @param fold_threshold numeric, default 2
#'   Minimum symmetric fold-difference required for a metabolite to be returned.
#' @param na_rm logical, default TRUE
#'   Whether to remove missing values when computing means.
#' @param pseudocount numeric, default 1e-8
#'   Small value added to means before division to avoid divide-by-zero.
#'
#' @return character
#'   Vector of metabolite column names meeting the fold-difference threshold.
#'
#' @examples
#' flagged_metabs <- get_metabs_2fold_vs_qc(df)
#'
#' @noRd
get_metabs_2fold_vs_qc <- function(
    df,
    metab_cols = NULL,
    class_col = "class",
    qc_label = "QC",
    fold_threshold = 2,
    na_rm = TRUE,
    pseudocount = 1e-8
) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  
  if (!class_col %in% names(df)) {
    stop(sprintf("`df` must contain a '%s' column.", class_col))
  }
  
  if (!is.numeric(fold_threshold) || length(fold_threshold) != 1L || fold_threshold < 1) {
    stop("`fold_threshold` must be a single numeric value >= 1.")
  }
  
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L || pseudocount < 0) {
    stop("`pseudocount` must be a single numeric value >= 0.")
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
  
  is_qc <- df[[class_col]] == qc_label
  is_sample <- !is.na(df[[class_col]]) & df[[class_col]] != qc_label
  
  if (!any(is_qc)) {
    stop(sprintf("No rows found where `%s == \"%s\"`.", class_col, qc_label))
  }
  
  if (!any(is_sample)) {
    stop(sprintf("No sample rows found where `%s != \"%s\"`.", class_col, qc_label))
  }
  
  sample_df <- df[is_sample, metab_cols, drop = FALSE]
  qc_df <- df[is_qc, metab_cols, drop = FALSE]
  
  non_numeric <- metab_cols[!vapply(sample_df, is.numeric, logical(1L))]
  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "All metabolite columns must be numeric. Non-numeric columns: %s",
        paste(non_numeric, collapse = ", ")
      )
    )
  }
  
  sample_means <- colMeans(sample_df, na.rm = na_rm)
  qc_means <- colMeans(qc_df, na.rm = na_rm)
  
  sample_means_adj <- sample_means + pseudocount
  qc_means_adj <- qc_means + pseudocount
  
  fold_diff <- pmax(
    sample_means_adj / qc_means_adj,
    qc_means_adj / sample_means_adj
  )
  
  names(fold_diff)[fold_diff >= fold_threshold]
}