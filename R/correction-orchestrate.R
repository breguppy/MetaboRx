#' Orchestrates correction step
#'
#' @keywords internal
#' @noRd
correct_data <- function(df, metab_cols, corMethod) {
  if (corMethod == "RF") {
    if (!requireNamespace("randomForest", quietly = TRUE))
      stop("Install 'randomForest' to use RF correction.", call. = FALSE)
    correction_str <- "Random Forest"
    parameters <- "builds 3 models using seeds 42, 31416, 272. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      return (rf_correction(df, metab_cols, ntree = 500, seed = seed))
      
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- .median_across_models(df_list, metadata_cols)
  } else if (corMethod == "LOESS") {
    correction_str <- "local polynomial regression"
    parameters <- "builds local polynomials of degree 2 that span 0.75 of the total QC values."
    df_corrected <- loess_correction(df, metab_cols, degree = 2)
  } else if (corMethod == "LC") {
    correction_str <- "local constant regression"
    parameters <- "fits study samples to the weighted mean of the nearby QC points."
    df_corrected <- loess_correction(df, metab_cols, degree = 0)
  } else if (corMethod == "LL") {
    correction_str <- "local linear regression"
    parameters <- "fits study samples to the line created by nearby QC points."
    df_corrected <- loess_correction(df, metab_cols, degree = 1)
  }else if (corMethod == "BW_RF") {
    if (!requireNamespace("randomForest", quietly = TRUE))
      stop("Install 'randomForest' to use RF correction.", call. = FALSE)
    correction_str <- "Batchwise Random Forest"
    parameters <- "3 models using seeds 42, 31416, 272 are built for each metabolite for each batch. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      return (bw_rf_correction(df, metab_cols, ntree = 500, seed = seed))
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- .median_across_models(df_list, metadata_cols)
  } else if (corMethod == "BW_LOESS") {
    correction_str <- "Batchwise LOESS"
    parameters <- "builds local polynomials of degree 2 that span 0.75 of the total QC values in each batch."
    df_corrected <- bw_loess_correction(df, metab_cols, degree = 2)
  }
  
  return(list(
    df = df_corrected,
    str = correction_str,
    parameters =  parameters
  ))
}