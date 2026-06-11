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

  tryCatch(
    {
      if (cor_method %in% c("Random Forest", "Batchwise Random Forest")) {
        met_scatter_rf(df_raw, df_cor, i = met_col)
      } else if (cor_method %in% c("local polynomial regression", "Batchwise LOESS", "local constant regression", "local linear regression")) {
        met_scatter_loess(df_raw, df_cor, cor_method, i = met_col)
      } else {
        ggplot2::ggplot() +
          ggplot2::labs(title = "No correction method selected.")
      }
    },
    error = function(e) {
      shiny::showNotification(paste("Scatter failed:", e$message),
        type = "error",
        duration = 8
      )
      ggplot2::ggplot() +
        ggplot2::labs(title = "Scatter failed \u2013 see notification")
    }
  )
}

#' @keywords internal
#' @noRd
# Helper for creating the RSD plot
make_rsd_plot <- function(p, d) {
  p <- .normalize_rsd_plot_params(p)

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
    shiny::need(ncol(df_after) > 4L, "No metabolites left after correction.")
  )

  tryCatch(
    {
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
    },
    error = function(e) {
      shiny::showNotification(
        paste("RSD comparison failed:", e$message),
        type = "error",
        duration = 8
      )
      ggplot2::ggplot() +
        ggplot2::labs(title = "RSD comparison failed \u2013 see notification")
    }
  )
}

#' @keywords internal
#' @noRd
.normalize_rsd_choice <- function(value,
                                  arg,
                                  choices,
                                  aliases = character()) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(
      "`", arg, "` must be one of: ",
      paste(choices, collapse = ", "),
      call. = FALSE
    )
  }

  if (value %in% names(aliases)) {
    value <- unname(aliases[[value]])
  }

  if (!value %in% choices) {
    stop(
      "`", arg, "` must be one of: ",
      paste(choices, collapse = ", "),
      call. = FALSE
    )
  }

  value
}

#' @keywords internal
#' @noRd
.normalize_rsd_plot_params <- function(p) {
  p$rsd_compare <- .normalize_rsd_choice(
    p$rsd_compare,
    "rsd_compare",
    c("filtered_cor_data", "transformed_cor_data")
  )
  p$rsd_plot_type <- .normalize_rsd_choice(
    p$rsd_plot_type,
    "rsd_plot_type",
    c("dist", "scatter")
  )
  p$rsd_cal <- .normalize_rsd_choice(
    p$rsd_cal,
    "rsd_cal",
    c("class_met", "met"),
    aliases = c("class-met" = "class_met")
  )

  p
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
    rsd_cal = c("class_met", "met"),
    rsd_compare = "filtered_cor_data",
    stringsAsFactors = FALSE
  )

  # Add transformed_cor_data variants only when transform is not "none"
  if (!identical(p$transform, "none")) {
    specs_trans <- specs
    specs_trans$rsd_compare <- "transformed_cor_data"
    specs <- rbind(specs, specs_trans)
  }

  rsd_plots <- vector("list", nrow(specs))
  plot_names <- character(nrow(specs))

  for (i in seq_len(nrow(specs))) {
    temp_params <- p
    temp_params$rsd_plot_type <- specs$rsd_plot_type[i]
    temp_params$rsd_compare <- specs$rsd_compare[i]
    temp_params$rsd_cal <- specs$rsd_cal[i]

    rsd_plots[[i]] <- make_rsd_plot(temp_params, d)
    plot_names[i] <- build_name(
      temp_params$rsd_plot_type,
      temp_params$rsd_compare,
      temp_params$rsd_cal
    )
  }

  list(
    rsd_plots   = rsd_plots,
    plot_names  = plot_names
  )
}
