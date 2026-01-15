#' Plot entry points used by the app
#' @keywords internal
#' @noRd
# Helper for creating metabolite scatter plot
make_met_scatter <- function(d, p, met_col) {
  # choose the correct plotting function based on the correction method.
  cor_method <- d$corrected$str
  df_raw <- d$filtered$df
  if (isTRUE(p$remove_imputed)) {
    df_cor <- d$filtered_corrected$df_mv
  } else {
    df_cor <- d$filtered_corrected$df_no_mv
  }
  
  tryCatch({
    if (cor_method %in% c("Random Forest", "Batchwise Random Forest")) {
      met_scatter_rf(df_raw, df_cor, i = met_col)
    } else if (cor_method %in% c("LOESS", "Batchwise LOESS")) {
      met_scatter_loess(df_raw, df_cor, i = met_col)
    } else {
      ggplot2::ggplot() + ggplot2::labs(title = "No correction method selected.")
    }
  }, error = function(e) {
    shiny::showNotification(paste("Scatter failed:", e$message),
                            type = "error",
                            duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "Scatter failed \u2013 see notification")
  })
}

#' @keywords internal
#' @noRd
# Helper for creating the RSD plot
make_rsd_plot <- function(p, d) {
  df_before <- d$filtered$df
  # Determine df_after based on rsd_compare selected by user.
  if (p$rsd_compare == "filtered_cor_data") {
    compared_to <- "Correction"
    if (isTRUE(p$remove_imputed)) {
      df_after <- d$filtered_corrected$df_mv
    } else {
      df_after <- d$filtered_corrected$df_no_mv
    }
  } else {
    compared_to <- "Correction and Transformation"
    if (isTRUE(p$remove_imputed)) {
      df_after <- d$transformed$df_mv
    } else {
      df_after <- d$transformed$df_no_mv
    }
  }
  
  # Need at least 1 metabolite column
  shiny::validate(
    shiny::need(ncol(df_before) > 4L, "No metabolites left before correction."),
    shiny::need(ncol(df_after)  > 4L, "No metabolites left after correction.")
  )
  
  tryCatch({
    if (identical(p$rsd_plot_type, "scatter")) {
      if (identical(p$rsd_cal, "met")) {
        plot_rsd_comparison(df_before, df_after, compared_to)
      } else {
        plot_rsd_comparison_class_met(df_before, df_after, compared_to)
      }
    } else {
      if (identical(p$rsd_cal, "met")) {
        plot_met_rsd_distributions(df_before, df_after, compared_to)
      } else {
        plot_class_rsd_distributions(df_before, df_after, compared_to)
      }
    }
  }, error = function(e) {
    shiny::showNotification(
      paste("RSD comparison failed:", e$message),
      type = "error",
      duration = 8
    )
    ggplot2::ggplot() + ggplot2::labs(title = "RSD comparison failed \u2013 see notification")
  })
}

#' @keywords internal
#' @noRd
# Helper for creating the PCA plot
make_pca_plot <- function(p, d) {
  # get after based on pca_compare selected by user.
  if (p$pca_compare == "filtered_cor_data") {
    after <- d$filtered_corrected$df_no_mv
    compared_to <- "Correction"
  } else {
    after <- d$transformed$df_no_mv
    compared_to <- "Correction and Transformation"
  }
  mets <- setdiff(names(after), c("sample", "batch", "class", "order"))
  shiny::validate(
    shiny::need(length(mets) >= 2, "Need at least 2 metabolite columns for PCA."),
    shiny::need(nrow(after) >= 3, "Need at least 3 samples for PCA.")
  )
  
  # Also ensure non-constant / non-NA columns
  X <- after[, mets, drop = FALSE]
  keep <- vapply(X, function(v) {
    v <- suppressWarnings(as.numeric(v))
    ok <- all(is.finite(v))
    nz <- (length(unique(v)) >= 2)
    ok && nz
  }, logical(1))
  shiny::validate(need(
    any(keep),
    "All metabolite columns are constant/invalid after filtering."
  ))
  
  # before data cannot have any missing values.
  before <- d$imputed$df
  tryCatch({
    plot_pca(p, before, after, compared_to)
  }, error = function(e) {
    shiny::showNotification(paste("PCA failed:", e$message),
                            type = "error",
                            duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "PCA failed \u2013 see notification")
  })
}

#' @keywords internal
#' @noRd
# Helper for creating the PCA plot
make_pca_loading_plot <- function(p, d) {
  # get after based on pca_compare selected by user.
  if (p$pca_compare == "filtered_cor_data") {
    after <- d$filtered_corrected$df_no_mv
    compared_to <- "Correction"
  } else {
    after <- d$transformed$df_no_mv
    compared_to <- "Correction and Transformation"
  }
  mets <- setdiff(names(after), c("sample", "batch", "class", "order"))
  shiny::validate(
    shiny::need(length(mets) >= 2, "Need at least 2 metabolite columns for PCA."),
    shiny::need(nrow(after) >= 3, "Need at least 3 samples for PCA.")
  )
  
  # Also ensure non-constant / non-NA columns
  X <- after[, mets, drop = FALSE]
  keep <- vapply(X, function(v) {
    v <- suppressWarnings(as.numeric(v))
    ok <- all(is.finite(v))
    nz <- (length(unique(v)) >= 2)
    ok && nz
  }, logical(1))
  shiny::validate(need(
    any(keep),
    "All metabolite columns are constant/invalid after filtering."
  ))
  
  # before data cannot have any missing values.
  before <- d$imputed$df
  tryCatch({
    plot_pca_loading(p, before, after, compared_to)
  }, error = function(e) {
    shiny::showNotification(paste("PCA failed:", e$message),
                            type = "error",
                            duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "PCA Loading failed \u2013 see notification")
  })
}

#' Make Hotelling PCA plot for report
#'
#' @param p List of parameters (must include qcImputeM, samImputeM, remove_imputed).
#' @param d List of data reactives/artifacts (must include filtered_corrected).
#'
#' @return A ggplot object or NULL if not available.
#'
#' @keywords internal
#' @noRd
make_hotelling_pca_plot <- function(p, d) {
  .require_pkg("ggplot2", "Hotelling PCA plot")
  
  if (is.null(d$filtered_corrected)) return(NULL)
  
  # Choose df version consistent with the app behavior
  remove_imputed <- isTRUE(p$remove_imputed)
  
  # In your server you call detect_hotelling_nonqc_dual_z() on df_no_mv
  # (even when remove_imputed is TRUE). Keep the same behavior here unless
  # you intentionally want it to differ.
  df <- d$filtered_corrected$df_no_mv
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  
  # Ensure p has the fields your detect_hotelling_nonqc_dual_z() expects
  p_hot <- list(
    qcImputeM  = p$qcImputeM %||% "nothing_to_impute",
    samImputeM = p$samImputeM %||% "nothing_to_impute"
  )
  
  res <- detect_hotelling_nonqc_dual_z(df, p_hot)
  
  if (!is.null(res$pca_plot)) res$pca_plot else NULL
}
