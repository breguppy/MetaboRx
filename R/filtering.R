#' Metabolite/value filtering functions
#'
#' @keywords internal
#' @noRd
filter_by_missing <- function(df, metab_cols, mv_cutoff) {
  # get metadata columns
  meta_cols <- setdiff(names(df), metab_cols)
  
  # Compute percentage of missing values per metabolite column
  missing_pct <- sapply(df[metab_cols], function(col) {
    mean(is.na(col)) * 100
  })
  
  # Keep only columns with missing percentage <= mv_cutoff
  mv_keep_cols <- metab_cols[missing_pct <= mv_cutoff]
  
  # list columns removed due to missing value %
  mv_removed_cols <- setdiff(metab_cols, mv_keep_cols)
  
  # Get all class-metabolite pairs where all values are missing
  classes_seen <- sort(unique(df$class))
  
  class_metab_all_missing <- if (length(classes_seen) == 0L || length(mv_keep_cols) == 0L) {
    data.frame(class = character(0), metabolite = character(0), n_rows_in_class = integer(0))
  } else {
    out <- vector("list", length(classes_seen))
    names(out) <- classes_seen
    
    for (cl in classes_seen) {
      idx <- which(df$class == cl)
      n_in_class <- length(idx)
      
      # If no rows for this class, treat as "all missing" for every kept metabolite
      if (n_in_class == 0L) {
        miss_all <- rep(TRUE, length(mv_keep_cols))
      } else {
        sub <- df[idx, mv_keep_cols, drop = FALSE]
        miss_all <- vapply(sub, function(x) all(is.na(x)), logical(1L))
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
      data.frame(class = character(0), metabolite = character(0), n_rows_in_class = integer(0))
  }
  
  # filter data by metabolite missing value
  df_filtered <- df[, c(meta_cols, mv_keep_cols), drop = FALSE]
  
  # Get metabolites that have QCs with missing values
  qc_idx <- which(df_filtered$class == "QC")
  if (length(mv_keep_cols) == 0L || length(qc_idx) == 0L) {
    qc_missing_mets <- character(0)
  } else {
    sub <- df_filtered[qc_idx, mv_keep_cols, drop = FALSE]
    qc_missing_mets <- mv_keep_cols[colSums(is.na(sub)) > 0]
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