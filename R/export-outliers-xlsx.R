#' Export outliers Excel file
#' @keywords internal
#' @noRd
export_outliers_xlsx <- function(p, d, file = NULL, hotelling_res = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  
  res <- hotelling_res
  if (is.null(res)) {
    df <- d$filtered_corrected$df_no_mv
    res <- detect_hotelling_nonqc_dual_z(df, p)
  }
  
  qc_label <- res$params$qc_label
  if (is.null(qc_label)) {
    qc_label <- "QC"
  }
  
  class_col <- res$params$class_col
  if (is.null(class_col)) {
    class_col <- "class"
  }
  
  value_z_threshold <- res$params$value_z_threshold
  flag_fraction <- res$params$flag_fraction
  
  if (is.null(value_z_threshold) || is.null(flag_fraction)) {
    stop(
      "The detector result must contain value_z_threshold and flag_fraction ",
      "in params. Recompute it using the updated detection function."
    )
  }
  
  ev <- res$extreme_values
  required_cols <- c("sample", class_col, "metabolite", "z_global", "T2")
  missing_cols <- setdiff(required_cols, names(ev))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in extreme_values: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  abs_z <- abs(ev$z_global)
  idx <- order(-abs_z, -ev$T2)
  
  candidate_values <- data.frame(
    Sample = ev$sample[idx],
    Class = ev[[class_col]][idx],
    Metabolite = ev$metabolite[idx],
    Z_score = ev$z_global[idx],
    Absolute_Z_score = abs_z[idx],
    Squared_Mahalanobis_distance = ev$T2[idx],
    stringsAsFactors = FALSE
  )
  
  study_values <- candidate_values[
    which(candidate_values$Class != qc_label),
    ,
    drop = FALSE
  ]
  qc_values <- candidate_values[
    which(candidate_values$Class == qc_label),
    ,
    drop = FALSE
  ]
  
  outlier_rows <- which(res$data$is_outlier_sample)
  outlier_samples <- data.frame(
    Sample = res$data$sample[outlier_rows],
    Class = res$data[[class_col]][outlier_rows],
    Squared_Mahalanobis_distance = res$data$T2[outlier_rows],
    stringsAsFactors = FALSE
  )
  outlier_samples <- outlier_samples[
    order(-outlier_samples$Squared_Mahalanobis_distance),
    ,
    drop = FALSE
  ]
  
  transformation_text <- if (isTRUE(res$params$log_transform)) {
    sprintf(
      "Z-scores use log2(intensity + %g).",
      res$params$log_offset
    )
  } else {
    "Z-scores use untransformed intensities."
  }
  
  selection_text <- sprintf(
    paste(
      "Values must first have an absolute Z-score greater than %g.",
      "The largest %g%% of qualifying absolute Z-scores across all metabolites",
      "in this sample pool are reported."
    ),
    value_z_threshold,
    100 * flag_fraction
  )
  
  selection_details <- paste(
    "The number reported is rounded down.",
    "Ties at the selection boundary are resolved by sample row,",
    "then metabolite column in the input data.",
    "Selection is independent of the sample's Mahalanobis flag."
  )
  
  wb <- openxlsx::createWorkbook()
  styles <- .xlsx_export_styles()
  
  shiny::withProgress(
    message = sprintf("Creating extreme_values_%s.xlsx...", Sys.Date()),
    value = 0,
    {
      .xlsx_write_described_sheet(
        wb = wb,
        sheet_name = "Potential Outlier Samples",
        description = paste(
          "Tab 1. Study and QC samples whose squared Mahalanobis distance",
          "exceeds the chi-square cutoff in PC1-PC2 space.",
          "PCA is fitted on study samples, and the ellipse is estimated",
          "from complete projected study and QC samples."
        ),
        x = outlier_samples,
        merge_cols = 12L,
        styles = styles
      )
      shiny::incProgress(
        1 / 3,
        detail = "Saved: Potential Outlier Samples"
      )
      
      .xlsx_write_described_sheet(
        wb = wb,
        sheet_name = "Study Extreme Values",
        description = paste(
          "Tab 2. Selected candidate extreme values in study samples.",
          "For each metabolite, Z-scores are calculated using the mean",
          "and sample standard deviation of all non-QC samples pooled together.",
          transformation_text,
          "Only finite observed values contribute to the Z-scores.",
          selection_text,
          selection_details
        ),
        x = study_values,
        merge_cols = 12L,
        styles = styles
      )
      shiny::incProgress(
        1 / 3,
        detail = "Saved: Study Extreme Values"
      )
      
      .xlsx_write_described_sheet(
        wb = wb,
        sheet_name = "QC Extreme Values",
        description = paste(
          "Tab 3. Selected candidate extreme values in QC samples.",
          "For each metabolite, Z-scores are calculated using the mean",
          "and sample standard deviation of QC samples only.",
          transformation_text,
          "Only finite observed values contribute to the Z-scores.",
          selection_text,
          selection_details
        ),
        x = qc_values,
        merge_cols = 12L,
        styles = styles
      )
      shiny::incProgress(
        1 / 3,
        detail = "Saved: QC Extreme Values"
      )
    }
  )
  
  .xlsx_save_or_return(wb, file)
}