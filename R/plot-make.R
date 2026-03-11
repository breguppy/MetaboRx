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
    } else if (cor_method %in% c("local polynomial regression", "Batchwise LOESS", "local constant regression", "local linear regression")) {
      met_scatter_loess(df_raw, df_cor, cor_method, i = met_col)
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
make_all_rsd_plots <- function(p, d) {
  build_name <- function(plot_type, compare, cal) {
    sprintf("rsd_%s_%s_%s", plot_type, compare, cal)
  }
  
  # Base grid: plot types x calcs (always for filtered_cor_data)
  specs <- expand.grid(
    rsd_plot_type = c("dist", "scatter"),
    rsd_cal       = c("class-met", "met"),
    rsd_compare   = "filtered_cor_data",
    stringsAsFactors = FALSE
  )
  
  # Add transformed_cor_data variants only when transform is not "none"
  if (!identical(p$transform, "none")) {
    specs_trans <- specs
    specs_trans$rsd_compare <- "transformed_cor_data"
    specs <- rbind(specs, specs_trans)
  }
  
  rsd_plots   <- vector("list", nrow(specs))
  plot_names  <- character(nrow(specs))
  
  for (i in seq_len(nrow(specs))) {
    temp_params <- p
    temp_params$rsd_plot_type <- specs$rsd_plot_type[i]
    temp_params$rsd_compare   <- specs$rsd_compare[i]
    temp_params$rsd_cal       <- specs$rsd_cal[i]
    
    rsd_plots[[i]]  <- make_rsd_plot(temp_params, d)
    plot_names[i]   <- build_name(temp_params$rsd_plot_type,
                                  temp_params$rsd_compare,
                                  temp_params$rsd_cal)
  }
  
  list(
    rsd_plots   = rsd_plots,
    plot_names  = plot_names
  )
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
