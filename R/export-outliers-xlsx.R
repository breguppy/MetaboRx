#' Export outliers Excel file
#'
#' @param p MetaboRx preprocessing parameters.
#' @param d MetaboRx data object.
#' @param file Optional destination for the workbook.
#' @param hotelling_res Optional cached result from `detect_outliers()`.
#'
#' @keywords internal
#' @noRd
export_outliers_xlsx <- function(p, d, file = NULL, hotelling_res = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  
  res <- hotelling_res
  if (is.null(res)) {
    df <- d$filtered_corrected$df_no_mv
    res <- detect_outliers(df, p)
  }
  
  params <- res$params
  required_params <- c(
    "min_rout_observations",
    "grubbs_alpha",
    "rout_q",
    "log_offset"
  )
  missing_params <- required_params[
    vapply(required_params, function(name) {
      is.null(params[[name]])
    }, logical(1))
  ]
  
  if (length(missing_params) > 0L) {
    stop(
      "Missing detector parameters: ",
      paste(missing_params, collapse = ", "),
      ". Recompute the result using the updated detect_outliers()."
    )
  }
  
  qc_label <- params$qc_label
  if (is.null(qc_label)) {
    qc_label <- "QC"
  }
  
  class_col <- params$class_col
  if (is.null(class_col)) {
    class_col <- "class"
  }
  
  ev <- res$extreme_values
  required_cols <- c(
    "sample",
    class_col,
    "metabolite",
    "method",
    "test_statistic",
    "adjusted_p",
    "p_cutoff",
    "T2"
  )
  missing_cols <- setdiff(required_cols, names(ev))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in extreme_values: ",
      paste(missing_cols, collapse = ", "),
      ". Recompute the result using the updated detect_outliers()."
    )
  }
  
  ev <- ev[order(ev$adjusted_p, -ev$T2), , drop = FALSE]
  
  candidate_values <- data.frame(
    Sample = ev$sample,
    Class = ev[[class_col]],
    Metabolite = ev$metabolite,
    Method = ev$method,
    Test_statistic = ev$test_statistic,
    Adjusted_p_value = ev$adjusted_p,
    P_value_cutoff = ev$p_cutoff,
    Squared_Mahalanobis_distance = ev$T2,
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
  
  method_description <- sprintf(
    paste(
      "For each metabolite, classes with fewer than %d observed values",
      "were screened, where estimable, using a two-sided Grubbs-type",
      "test (alpha = %g). The test used externally studentized one-way",
      "ANOVA residuals with within-class variance pooled across all classes",
      "for that metabolite and p-values Bonferroni-adjusted within each class.",
      "Classes with %d or more observed values were screened using",
      "ROUT (Q = %g%%). Both tests used finite observed values transformed",
      "as log2(intensity + %g). Flagged values were retained.",
      "Value screening was independent of the sample's Mahalanobis flag."
    ),
    params$min_rout_observations,
    params$grubbs_alpha,
    params$min_rout_observations,
    100 * params$rout_q,
    params$log_offset
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
          "Tab 2. Candidate extreme values in study samples.",
          method_description
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
          "Tab 3. Candidate extreme values in QC samples.",
          method_description
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