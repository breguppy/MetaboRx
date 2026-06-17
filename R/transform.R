#' Internal standard normalization (row-wise)
#'
#' For each row, compute the mean of internal standard columns (ISTD*/ITSD*),
#' then divide all non-ISTD metabolite columns by that mean.
#'
#' @param df            Data frame containing metadata + metabolite columns.
#' @param metab_cols    Character vector of metabolite columns to normalize (non-ISTD).
#' @param istd_cols     Character vector of internal standard columns.
#' @param min_istd      Minimum number of non-missing ISTD values required to compute a mean.
#' @param na_action     What to do if ISTD mean cannot be computed for a row.
#'                      One of c("leave", "na", "error").
#'
#' @return `df` with `metab_cols` normalized in-place.
#'
#' @keywords internal
#' @noRd
.istd_norm <- function(df,
                       metab_cols,
                       istd_cols,
                       min_istd = 1L,
                       na_action = c("leave", "na", "error")) {
  na_action <- match.arg(na_action)

  if (length(istd_cols) == 0L) {
    stop("ISTD_norm requested but no ISTD/ITSD columns were found.")
  }
  if (length(metab_cols) == 0L) {
    return(df)
  }

  istd_data <- df[, istd_cols, drop = FALSE]

  # Row-wise ISTD mean, requiring at least `min_istd` non-missing values.
  n_nonmiss <- rowSums(!is.na(istd_data))
  istd_mean <- rowMeans(istd_data, na.rm = TRUE)
  istd_mean[n_nonmiss < min_istd] <- NA_real_

  if (na_action == "error" && anyNA(istd_mean)) {
    bad_rows <- which(is.na(istd_mean))
    stop(sprintf(
      "Cannot compute ISTD mean for %d row(s) (need >= %d non-missing ISTD values). Example row(s): %s",
      length(bad_rows),
      min_istd,
      paste(utils::head(bad_rows, 10), collapse = ", ")
    ))
  }

  metab_data <- df[, metab_cols, drop = FALSE]

  # Divide each row by its ISTD mean
  norm_data <- sweep(metab_data, 1, istd_mean, FUN = "/")

  if (na_action == "leave") {
    # For rows where ISTD mean is NA, keep original values
    bad <- is.na(istd_mean)
    if (any(bad)) norm_data[bad, ] <- metab_data[bad, ]
  } else if (na_action == "na") {
    # For rows where ISTD mean is NA, force normalized values to NA
    bad <- is.na(istd_mean)
    if (any(bad)) norm_data[bad, ] <- NA_real_
  }

  df[, metab_cols] <- norm_data
  df
}

#' Transformation methods for corrected data
#'
#' @keywords internal
#' @noRd
.total_ratio_norm <- function(df, metab_cols) {
  metab_data <- df[, metab_cols, drop = FALSE]

  # sum metab_cols values in each row (sample).
  row_sums <- rowSums(metab_data, na.rm = TRUE)
  # compute number of non-missing values in each row (sample).
  non_missing_counts <- rowSums(!is.na(metab_data))

  # determine row ratio = (number of non-missing) / (row sum)
  ratios <- non_missing_counts / row_sums

  # multiply each row by ratio.
  trn_data <- sweep(metab_data, 1, ratios, FUN = "*")

  df[, metab_cols] <- trn_data

  return(df)
}


transform_data <- function(filtered_corrected, transform, withheld_cols, ex_ISTD = TRUE) {
  df_mv <- filtered_corrected$df_mv
  df_no_mv <- filtered_corrected$df_no_mv

  meta_cols <- c("sample", "batch", "class", "order")
  metab_cols_mv <- setdiff(names(df_mv), meta_cols)
  metab_cols_no_mv <- setdiff(names(df_no_mv), meta_cols)

  # Identify ISTD columns from the full metabolite sets
  istd_names_mv <- metab_cols_mv[grepl("ISTD|ITSD", metab_cols_mv, ignore.case = FALSE)]
  istd_names_no_mv <- metab_cols_no_mv[grepl("ISTD|ITSD", metab_cols_no_mv, ignore.case = FALSE)]

  withheld_cols_mv <- withheld_cols
  withheld_cols_no_mv <- withheld_cols

  # Exclude ISTDs from TRN calculation when requested
  if (isTRUE(ex_ISTD)) {
    withheld_cols_mv <- unique(c(withheld_cols_mv, istd_names_mv))
    withheld_cols_no_mv <- unique(c(withheld_cols_no_mv, istd_names_no_mv))
  }

  # Columns to transform
  metab_cols_mv_trn <- setdiff(metab_cols_mv, withheld_cols_mv)
  metab_cols_no_mv_trn <- setdiff(metab_cols_no_mv, withheld_cols_no_mv)

  transformed_df_mv <- df_mv
  transformed_df_no_mv <- df_no_mv

  equal_weight_df_mv <- transformed_df_mv
  equal_weight_df_no_mv <- transformed_df_no_mv

  if (transform == "none") {
    transform_str <- "After correction, no scaling or transformations have been applied."
  } else if (transform == "ISTD_norm") {
    transform_str <- paste(
      "After correction, all metabolite levels are normalized to the average of",
      "the internal standards within that sample."
    )

    metab_cols_mv_norm <- setdiff(metab_cols_mv, withheld_cols_mv)
    metab_cols_no_mv_norm <- setdiff(metab_cols_no_mv, withheld_cols_no_mv)

    # Never normalize ISTD columns themselves
    metab_cols_mv_norm <- setdiff(metab_cols_mv_norm, istd_names_mv)
    metab_cols_no_mv_norm <- setdiff(metab_cols_no_mv_norm, istd_names_no_mv)

    transformed_df_mv <- .istd_norm(
      transformed_df_mv,
      metab_cols = metab_cols_mv_norm,
      istd_cols = istd_names_mv,
      min_istd = 1L,
      na_action = "leave"
    )

    transformed_df_no_mv <- .istd_norm(
      transformed_df_no_mv,
      metab_cols = metab_cols_no_mv_norm,
      istd_cols = istd_names_no_mv,
      min_istd = 1L,
      na_action = "leave"
    )
  } else if (transform == "TRN") {
    transform_str <- paste(
      "After correction, metabolite level values are ratiometrically normalized",
      "to total metabolite signal on a per sample basis. This normalization is",
      "done by summing all individual post-QC corrected metabolite level values",
      "within a sample (total signal) and then dividing each individual",
      "metabolite level value within that sample by the total signal. Next, the",
      "individual values are multiplied by the total number of metabolites",
      "present in the sample for easier visualization. This normalization",
      "quantifies individual metabolite values across samples based on their",
      "proportion to total metabolite load, in arbitrary units, within each",
      "individual sample. Because arbitrary units for a given metabolite",
      "quantitatively scale across samples, levels of a given metabolite may be",
      "quantitatively compared across samples. Because unit scaling is different",
      "for each metabolite, different metabolites within a sample cannot be",
      "quantitatively compared. However, because differences in arbitrary unit",
      "scaling between samples cancel out by division, within-sample metabolite",
      "ratios can be quantitatively compared across samples."
    )

    # Equal-weight only TRN-included metabolites using NON-QC samples
    equal_weight_df_mv <- equally_weight_metabolites(
      df = transformed_df_mv,
      metab_cols = metab_cols_mv_trn,
      class_col = "class",
      qc_label = "QC",
      target_mean = 1,
      na_rm = TRUE
    )

    equal_weight_df_no_mv <- equally_weight_metabolites(
      df = transformed_df_no_mv,
      metab_cols = metab_cols_no_mv_trn,
      class_col = "class",
      qc_label = "QC",
      target_mean = 1,
      na_rm = TRUE
    )

    # Apply TRN only to TRN-included metabolite columns
    transformed_df_mv <- .total_ratio_norm(equal_weight_df_mv, metab_cols_mv_trn)
    transformed_df_no_mv <- .total_ratio_norm(equal_weight_df_no_mv, metab_cols_no_mv_trn)

    # Remove ISTD columns from the returned equal-weight dfs
    equal_weight_df_mv <- equal_weight_df_mv[, setdiff(names(equal_weight_df_mv), istd_names_mv), drop = FALSE]
    equal_weight_df_no_mv <- equal_weight_df_no_mv[, setdiff(names(equal_weight_df_no_mv), istd_names_no_mv), drop = FALSE]
  } else if (transform == "PQN") {
    transform_str <- paste(
      "This tab shows normalized metabolite level values using probabilistic ",
      "quotient normalization (PQN). PQN is a sample-based normalization method ",
      "computed in 3 steps:(1) Each metabolite is divided by the median value ",
      "of that metabolite across all samples. (2) For each sample, the median ",
      "of these quotients is computed as an estimate of the sample's most ",
      "probable dilution factor. (3) Post-QC-corrected metabolite intensities ",
      "are divided by this sample-specific median quotient. This normalization ",
      "rescales each sample by its median fold difference relative to a reference ",
      "spectrum, correcting for global dilution or concentration differences ",
      "while preserving relative biological differences in individual ",
      "metabolites. Data remain in arbitrary units. Because arbitary units for ",
      "a given metabolite quantitatively scale across samples, levels of a given ",
      "metabolite may be quantiatively compared across samples. Because unit ",
      "scaling is different for each metabolite, different metabolites within ",
      "in a sample cannot be quantitatively compared. However, because ",
      "differences in arbitrary unit scaling between samples cancel out by ",
      "divsion, within-sample metabolite ratios can be quantitatively compared ",
      "across samples."
    )
    keep_cols_mv <- setdiff(names(transformed_df_mv), withheld_cols_mv)
    keep_cols_no_mv <- setdiff(names(transformed_df_no_mv), withheld_cols_no_mv)
    transformed_df_mv <- pqn_norm(
      df = transformed_df_mv[, keep_cols_mv],
      metab_cols = setdiff(metab_cols_mv, withheld_cols_mv)
    )
    transformed_df_no_mv <- pqn_norm(
      df = transformed_df_no_mv[, keep_cols_no_mv],
      metab_cols = setdiff(metab_cols_no_mv, withheld_cols_no_mv)
    )
  }

  return(list(
    df_mv = transformed_df_mv,
    df_no_mv = transformed_df_no_mv,
    equal_weight_df_mv = equal_weight_df_mv,
    equal_weight_df_no_mv = equal_weight_df_no_mv,
    str = transform_str,
    withheld_cols_mv = withheld_cols_mv,
    withheld_cols_no_mv = withheld_cols_no_mv
  ))
}

#' Equally weight metabolites using NON-QC samples
#'
#' Scales each metabolite column so that the mean across NON-QC samples
#' equals `target_mean` (default = 1000).
#'
#' @param df A data.frame containing a class column.
#' @param metab_cols Character vector of metabolite columns to scale.
#' @param class_col Name of class column.
#' @param qc_label Value in `class_col` identifying QC samples.
#' @param target_mean Desired mean across NON-QC samples.
#' @param na_rm Logical; whether to ignore NA values when computing means.
#'
#' @return A data.frame with scaled metabolite columns.
#'
#' @keywords internal
#' @noRd
equally_weight_metabolites <- function(df,
                                       metab_cols,
                                       class_col = "class",
                                       qc_label = "QC",
                                       target_mean = 1,
                                       na_rm = TRUE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }

  if (!class_col %in% names(df)) {
    stop(sprintf("`df` must contain a '%s' column.", class_col))
  }

  missing_cols <- setdiff(metab_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Columns not found in `df`: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (length(metab_cols) == 0L) {
    return(df)
  }

  non_qc_idx <- !is.na(df[[class_col]]) & df[[class_col]] != qc_label

  if (!any(non_qc_idx)) {
    stop("No non-QC samples found; cannot compute equal-weight scaling.")
  }

  metab_data <- as.data.frame(df[, metab_cols, drop = FALSE])

  # Force numeric once, up front
  metab_data[] <- lapply(metab_data, function(x) {
    if (is.numeric(x)) x else suppressWarnings(as.numeric(x))
  })

  # Means computed ONLY on non-QC rows
  mean_fun <- if (isTRUE(na_rm)) {
    function(x) mean(x, na.rm = TRUE)
  } else {
    function(x) mean(x, na.rm = FALSE)
  }

  col_means <- vapply(
    metab_data[non_qc_idx, , drop = FALSE],
    mean_fun,
    numeric(1L)
  )

  bad_cols <- names(col_means)[is.na(col_means) | abs(col_means) < .Machine$double.eps]

  scale_factors <- target_mean / col_means
  scale_factors[bad_cols] <- 1

  # Vectorized column-wise multiplication
  scaled_data <- sweep(
    as.matrix(metab_data),
    MARGIN = 2,
    STATS = scale_factors,
    FUN = "*"
  )

  df[, metab_cols] <- scaled_data

  return(df)
}
#' PQN normalization method
#'
#' Uses `pmp::pqn_normalisation()` on non-QC samples only.
#' `pmp::pqn_normalisation()` expects metabolites/features in rows and samples
#' in columns, so the metabolite matrix is transposed before normalization and
#' transposed back afterward.
#'
#' Scaling factors are computed as each sample's median ratio to the
#' metabolite-wise median across non-QC samples.
#'
#' @param df A data.frame containing metadata columns and metabolite columns.
#' @param metab_cols Character vector of metabolite column names.
#' @param class_col Name of the class column.
#' @param qc_label Label identifying QC samples.
#' @param na_rm Currently unused; included for API compatibility.
#'
#' @return A data.frame with the same columns and row order as `df`.
#'
#' @noRd
pqn_norm <- function(df,
                     metab_cols,
                     class_col = "class",
                     qc_label = "QC",
                     na_rm = TRUE) {
  if (!requireNamespace("pmp", quietly = TRUE)) {
    stop("Install 'pmp' to use PQN normalization.", call. = FALSE)
  }

  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.", call. = FALSE)
  }

  required_meta_cols <- c("sample", "batch", class_col, "order")
  missing_meta_cols <- setdiff(required_meta_cols, names(df))

  if (length(missing_meta_cols) > 0L) {
    stop(
      sprintf(
        "`df` must contain these metadata columns: %s",
        paste(missing_meta_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_metab_cols <- setdiff(metab_cols, names(df))

  if (length(missing_metab_cols) > 0L) {
    stop(
      sprintf(
        "Columns not found in `df`: %s",
        paste(missing_metab_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(metab_cols) == 0L) {
    return(df)
  }

  original_col_order <- names(df)
  original_row_order <- seq_len(nrow(df))

  work_df <- df
  work_df$.original_row_order <- original_row_order

  is_qc <- work_df[[class_col]] == qc_label
  is_qc[is.na(is_qc)] <- FALSE

  qc_df <- work_df[is_qc, , drop = FALSE]
  non_qc_df <- work_df[!is_qc, , drop = FALSE]

  if (nrow(non_qc_df) == 0L) {
    warning("No non-QC samples found. Returning `df` unchanged.", call. = FALSE)
    return(df)
  }

  if (anyDuplicated(non_qc_df$sample)) {
    sample_names <- make.unique(as.character(non_qc_df$sample))
  } else {
    sample_names <- as.character(non_qc_df$sample)
  }

  pqn_input <- t(as.matrix(non_qc_df[, metab_cols, drop = FALSE]))
  rownames(pqn_input) <- metab_cols
  colnames(pqn_input) <- sample_names

  classes <- rep("sample", ncol(pqn_input))

  suppressWarnings(
    pqn_data <- pmp::pqn_normalisation(
      df = pqn_input,
      classes = classes,
      qc_label = "all",
      ref_method = "median"
    )
  )

  pqn_data <- as.data.frame(t(pqn_data), check.names = FALSE)
  pqn_data$sample <- rownames(pqn_data)

  non_qc_meta <- non_qc_df[, setdiff(names(work_df), metab_cols), drop = FALSE]
  non_qc_meta$.pqn_sample_name <- sample_names

  normalized_non_qc_df <- merge(
    non_qc_meta,
    pqn_data,
    by.x = ".pqn_sample_name",
    by.y = "sample",
    all.x = TRUE,
    sort = FALSE
  )

  normalized_non_qc_df$.pqn_sample_name <- NULL

  combined_df <- rbind(
    normalized_non_qc_df[, names(work_df), drop = FALSE],
    qc_df[, names(work_df), drop = FALSE]
  )

  combined_df <- combined_df[order(combined_df$.original_row_order), , drop = FALSE]
  combined_df$.original_row_order <- NULL

  combined_df <- combined_df[, original_col_order, drop = FALSE]
  rownames(combined_df) <- NULL

  combined_df
}
