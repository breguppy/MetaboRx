#' Impute missing values with user selected options
#'
#' @keywords internal
#' @noRd
impute_missing <- function(df, metab_cols, qcImputeM, samImputeM) {
  # Count number of missing values to impute
  n_missv <- sum(is.na(df[, metab_cols]))
  imputed_df <- df

  # helper function to apply an imputation strategy to a subset of data
  apply_impute <- function(sub_df, method) {
    if (method == "nothing_to_impute") {
      sub_df <- sub_df
      str <- "nothing to impute"
    } else if (method == "median") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- median(col, na.rm = TRUE)
        col
      })
      str <- "metabolite median"
    } else if (method == "mean") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
      })
      str <- "metabolite mean"
    } else if (method == "class_median") {
      sub_df <- sub_df |>
        dplyr::group_by(.data[["class"]]) |>
        dplyr::mutate(dplyr::across(dplyr::all_of(metab_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) |>
        dplyr::ungroup()
      str <- "class-metabolite median"
    } else if (method == "class_mean") {
      sub_df <- sub_df |>
        dplyr::group_by(.data[["class"]]) |>
        dplyr::mutate(dplyr::across(dplyr::all_of(metab_cols), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) |>
        dplyr::ungroup()
      str <- "class-metabolite mean"
    } else if (method == "min") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- min(col, na.rm = TRUE)
        col
      })
      str <- "minimum metabolite value"
    } else if (method == "minHalf") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- 0.5 * min(col, na.rm = TRUE)
        col
      })
      str <- "half the minimum metabolite value"
    } else if (method == "zero") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- 0
        col
      })
      str <- "zero"
    } else if (method == "KNN") {
      metab_matrix <- as.matrix(sub_df[metab_cols])
      transposed <- t(metab_matrix)
      knn_result <- impute::impute.knn(transposed,
        rowmax = 0.99,
        colmax = 0.99,
        maxp = 15000
      )
      imputed_matrix <- t(knn_result$data)
      sub_df[metab_cols] <- as.data.frame(imputed_matrix)
      str <- "KNN"
    }
    return(list(sub_df = sub_df, str = str))
  }

  # Identify QC and Sample rows
  is_qc <- imputed_df$class == "QC"
  qc_df <- imputed_df[is_qc, ]
  sam_df <- imputed_df[!is_qc, ]

  # Apply user chosen imputation methods
  qc_imputed <- apply_impute(qc_df, qcImputeM)
  qc_df <- qc_imputed$sub_df
  qc_str <- qc_imputed$str

  sam_imputed <- apply_impute(sam_df, samImputeM)
  sam_df <- sam_imputed$sub_df
  sam_str <- sam_imputed$str

  # Combine the results and make sure its in order
  imputed_df <- dplyr::bind_rows(qc_df, sam_df) |>
    dplyr::arrange(.data[["order"]])

  return(list(
    df = imputed_df,
    qc_str = qc_str,
    sam_str = sam_str,
    n_missv = n_missv
  ))
}
