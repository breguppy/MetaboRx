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
  
  ev <- res$extreme_values
  required_cols <- c("sample", class_col, "metabolite", "z_class", "T2")
  missing_cols <- setdiff(required_cols, names(ev))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in extreme_values: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  abs_z_class <- abs(ev$z_class)
  idx <- order(-abs_z_class, -ev$T2)
  
  candidate_values <- data.frame(
    Sample = ev$sample[idx],
    Class = ev[[class_col]][idx],
    Metabolite = ev$metabolite[idx],
    z_class = ev$z_class[idx],
    abs_z_class = abs_z_class[idx],
    T2 = ev$T2[idx],
    stringsAsFactors = FALSE
  )
  
  study_values <- candidate_values[
    which(candidate_values$Class != qc_label), , drop = FALSE
  ]
  qc_values <- candidate_values[
    which(candidate_values$Class == qc_label), , drop = FALSE
  ]
  
  outlier_rows <- which(res$data$is_outlier_sample)
  outlier_samples <- data.frame(
    Sample = res$data$sample[outlier_rows],
    Class = res$data[[class_col]][outlier_rows],
    T2 = res$data$T2[outlier_rows],
    stringsAsFactors = FALSE
  )
  outlier_samples <- outlier_samples[
    order(-outlier_samples$T2), , drop = FALSE
  ]
  
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
      shiny::incProgress(1 / 3, detail = "Saved: Potential Outlier Samples")
      
      .xlsx_write_described_sheet(
        wb = wb,
        sheet_name = "Study Extreme Values",
        description = paste(
          "Tab 2. Candidate extreme metabolite values in study samples.",
          "A value is flagged when its absolute Z-score within its study",
          "class exceeds 3, regardless of the sample's ellipse status."
        ),
        x = study_values,
        merge_cols = 12L,
        styles = styles
      )
      shiny::incProgress(1 / 3, detail = "Saved: Study Extreme Values")
      
      .xlsx_write_described_sheet(
        wb = wb,
        sheet_name = "QC Extreme Values",
        description = paste(
          "Tab 3. Candidate extreme metabolite values in QC samples.",
          "A value is flagged when its absolute Z-score within the QC",
          "class exceeds 3, regardless of the sample's ellipse status."
        ),
        x = qc_values,
        merge_cols = 12L,
        styles = styles
      )
      shiny::incProgress(1 / 3, detail = "Saved: QC Extreme Values")
    }
  )
  
  .xlsx_save_or_return(wb, file)
}