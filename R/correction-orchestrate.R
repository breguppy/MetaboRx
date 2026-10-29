#' Orchestrates correction step
#'
#' @keywords internal
#' @noRd
correct_data <- function(df, metab_cols, corMethod, met_scaling) {
  if (corMethod %in% c("RF", "BW_RF")) {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Install 'randomForest' to use RF correction.", call. = FALSE)
    }
  }

  if (corMethod == "RF") {
    correction_str <- "Random Forest"
    parameters <- "builds 3 models using seeds 42, 31416, 272. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      rf_correction(df, metab_cols, ntree = 500, seed = seed, met_scaling = met_scaling)
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- .median_across_models(df_list, metadata_cols)
  } else if (corMethod == "LOESS") {
    correction_str <- "local polynomial regression"
    parameters <- "builds local polynomials for QC trends with adaptive safeguards; unstable degree-2 fits fall back to simpler local fits."
    df_corrected <- loess_correction(df, metab_cols, degree = 2, met_scaling = met_scaling)
  } else if (corMethod == "LC") {
    correction_str <- "local constant regression"
    parameters <- "fits study samples using a Nadaraya-Watson kernel-weighted mean of nearby QC points."
    df_corrected <- nw_correction(df, metab_cols, met_scaling = met_scaling, span = 0.75)
  } else if (corMethod == "LL") {
    correction_str <- "local linear regression"
    parameters <- "fits study samples to the line created by nearby QC points."
    df_corrected <- loess_correction(df, metab_cols, degree = 1, met_scaling = met_scaling)
  } else if (corMethod == "BW_RF") {
    correction_str <- "Batchwise Random Forest"
    parameters <- "3 models using seeds 42, 31416, 272 are built for each metabolite for each batch. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      bw_rf_correction(df, metab_cols, ntree = 500, seed = seed)
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- .median_across_models(df_list, metadata_cols)
  } else if (corMethod == "BW_LOESS") {
    correction_str <- "Batchwise LOESS"
    parameters <- "builds local polynomials for QC trends within each batch with adaptive safeguards; unstable degree-2 fits fall back to simpler local fits."
    df_corrected <- bw_loess_correction(df, metab_cols, degree = 2)
  } else if (corMethod == "None") {
    correction_str <- "none"
    parameters <- "does not apply any drift correction models."
    if (met_scaling == "original") {
      df_corrected <- df
    } else {
      qcid <- which(df$class == "QC")
      out <- df
      
      for (metab in metab_cols) {
        x <- df[[metab]]
        
        if (!is.numeric(x)) {
          stop(sprintf("Metabolite '%s' must be numeric.", metab))
        }
        
        qc_y <- x[qcid]
        valid_qc <- is.finite(qc_y) & qc_y > 0
        
        if (!any(valid_qc)) {
          warning(
            sprintf(
              "Skipping scaling for '%s': no finite positive QC values.",
              metab
            )
          )
          next
        }
        
        out[[metab]] <- x / stats::median(qc_y[valid_qc])
      }
      
      df_corrected <- out
    }
  }

  list(
    df = df_corrected,
    str = correction_str,
    parameters = parameters,
    met_scaling = met_scaling
  )
}
