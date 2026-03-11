# pca_utils.R
#
# Utility functions for computing, plotting, and exporting PCA results.
#
# This module centralizes PCA preparation so that:
#   - PCA score plots and loading plots reuse the same PCA fit
#   - PCA loading Excel exports use the same underlying PCA results
#   - repeated PCA computation is minimized during figure export
#
# Expected metadata columns are:
#   sample, batch, class, order
#
# @keywords internal

#' Get shared metabolite columns for paired PCA
#'
#' Identifies the overlapping non-metadata columns between two data frames.
#' These shared metabolite columns are used so that the "before" and "after"
#' PCA models are built from directly comparable feature sets.
#'
#' @param before data.frame
#'   First data frame.
#' @param after data.frame
#'   Second data frame.
#' @param meta_cols character
#'   Metadata columns to exclude.
#'
#' @return character
#'   Shared metabolite column names.
#'
#' @noRd
get_shared_pca_metab_cols <- function(
    before,
    after,
    meta_cols = c("sample", "batch", "class", "order")
) {
  intersect(
    setdiff(names(before), meta_cols),
    setdiff(names(after), meta_cols)
  )
}


#' Prepare PCA matrix
#'
#' Restricts the input data frame to a specified set of metabolite columns,
#' imputes missing values if needed, and validates that all PCA inputs are
#' numeric.
#'
#' @param df data.frame
#'   Input data frame containing metadata and metabolite columns.
#' @param p list
#'   Parameter list containing imputation settings.
#' @param metab_cols character
#'   Metabolite columns to include in PCA.
#'
#' @return data.frame
#'   Numeric data frame containing only metabolite columns.
#'
#' @noRd
prep_pca_matrix <- function(df, p, metab_cols) {
  if (length(metab_cols) == 0L) {
    stop("No metabolite columns available for PCA.")
  }
  
  out <- df[, metab_cols, drop = FALSE]
  
  if (anyNA(out)) {
    results <- impute_missing(df, metab_cols, p$qcImputeM, p$samImputeM)
    out <- results$df[, metab_cols, drop = FALSE]
  }
  
  is_num <- vapply(out, is.numeric, logical(1))
  if (!all(is_num)) {
    bad_cols <- names(out)[!is_num]
    stop(
      sprintf(
        "All PCA columns must be numeric. Non-numeric columns found: %s",
        paste(bad_cols, collapse = ", ")
      )
    )
  }
  
  out
}


#' Compute PCA for a single data frame
#'
#' Fits PCA once and returns the model, sample scores, variable loadings,
#' and explained variance.
#'
#' @param df data.frame
#'   Input data frame containing metadata and metabolite columns.
#' @param p list
#'   Parameter list containing imputation settings.
#' @param metab_cols character
#'   Metabolite columns to include in PCA.
#' @param meta_cols character
#'   Metadata columns to retain in the score output.
#'
#' @return list
#'   A list containing:
#'   \itemize{
#'     \item \code{fit}: prcomp result
#'     \item \code{scores}: sample scores with metadata appended
#'     \item \code{loadings}: variable loadings
#'     \item \code{explained_variance}: explained variance table
#'     \item \code{metab_cols}: metabolite columns used in the PCA
#'   }
#'
#' @noRd
compute_single_pca <- function(
    df,
    p,
    metab_cols,
    meta_cols = c("sample", "batch", "class", "order")
) {
  x <- prep_pca_matrix(df = df, p = p, metab_cols = metab_cols)
  fit <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  
  available_meta <- intersect(meta_cols, names(df))
  
  scores_df <- as.data.frame(fit$x, stringsAsFactors = FALSE)
  scores_df <- dplyr::bind_cols(scores_df, df[, available_meta, drop = FALSE])
  
  loadings_df <- as.data.frame(fit$rotation, stringsAsFactors = FALSE)
  loadings_df$variable <- rownames(loadings_df)
  rownames(loadings_df) <- NULL
  loadings_df <- loadings_df[, c("variable", setdiff(names(loadings_df), "variable")), drop = FALSE]
  
  explained_var <- (fit$sdev^2) / sum(fit$sdev^2)
  explained_variance_df <- data.frame(
    PC = paste0("PC", seq_along(explained_var)),
    explained_variance = as.numeric(explained_var),
    cumulative_explained_variance = cumsum(as.numeric(explained_var)),
    stringsAsFactors = FALSE
  )
  
  list(
    fit = fit,
    scores = scores_df,
    loadings = loadings_df,
    explained_variance = explained_variance_df,
    metab_cols = metab_cols
  )
}


#' Compute paired PCA results for before/after comparison
#'
#' Uses the same shared metabolite columns for both data sets so that the
#' resulting PCA summaries are directly comparable.
#'
#' @param before data.frame
#'   Data frame for the "Before" dataset.
#' @param after data.frame
#'   Data frame for the "After" dataset.
#' @param p list
#'   Parameter list containing imputation settings.
#' @param before_label character
#'   Label for the first dataset.
#' @param after_label character
#'   Label for the second dataset.
#' @param meta_cols character
#'   Metadata columns to retain in the score output.
#'
#' @return list
#'   Paired PCA results.
#'
#' @noRd
compute_pca_pair <- function(
    before,
    after,
    p,
    before_label = "Before",
    after_label = "After",
    meta_cols = c("sample", "batch", "class", "order")
) {
  metab_cols <- get_shared_pca_metab_cols(
    before = before,
    after = after,
    meta_cols = meta_cols
  )
  
  before_res <- compute_single_pca(
    df = before,
    p = p,
    metab_cols = metab_cols,
    meta_cols = meta_cols
  )
  
  after_res <- compute_single_pca(
    df = after,
    p = p,
    metab_cols = metab_cols,
    meta_cols = meta_cols
  )
  
  before_res$label <- before_label
  after_res$label <- after_label
  
  list(
    before = before_res,
    after = after_res,
    meta_cols = meta_cols,
    metab_cols = metab_cols
  )
}


#' Build PCA score plot from precomputed PCA results
#'
#' @param p list
#'   Parameter list.
#' @param pca_pair list
#'   Output from \code{compute_pca_pair()}.
#' @param compared_to character
#'   Text appended to the plot title.
#'
#' @return ggplot
#'   Combined PCA score plot.
#'
#' @noRd
plot_pca_from_result <- function(p, pca_pair, compared_to) {
  before_df <- pca_pair$before$scores
  after_df <- pca_pair$after$scores
  
  if (!all(c("PC1", "PC2") %in% names(before_df))) {
    stop("PCA scores do not contain PC1 and PC2 for the before dataset.")
  }
  if (!all(c("PC1", "PC2") %in% names(after_df))) {
    stop("PCA scores do not contain PC1 and PC2 for the after dataset.")
  }
  
  combined <- dplyr::bind_rows(before_df, after_df)
  x_limits <- range(combined$PC1, na.rm = TRUE)
  y_limits <- range(combined$PC2, na.rm = TRUE)
  
  var_raw <- 100 * pca_pair$before$explained_variance$explained_variance[1:2]
  var_cor <- 100 * pca_pair$after$explained_variance$explained_variance[1:2]
  
  cbPalette <- c(
    "#F3C300", "#875692", "#ee7733", "#A1CAF1", "#BE0032",
    "#C2B280", "#555555", "#008856", "#E68FAC", "#0067A5",
    "#F99379", "#332288", "#F6A600", "#B3446C", "#DCD300",
    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
    "#bbbbbb", "#000000", "#33bbee", "#ccddaa", "#225555"
  )
  
  col <- p$color_col %||% "class"
  use_gradient <- identical(col, "order")
  
  if (!col %in% names(before_df) || !col %in% names(after_df)) {
    stop(sprintf("Column '%s' not found in PCA score data.", col))
  }
  
  if (use_gradient) {
    combined[[col]] <- as.numeric(combined[[col]])
    before_df[[col]] <- as.numeric(before_df[[col]])
    after_df[[col]] <- as.numeric(after_df[[col]])
    
    order_range <- range(combined[[col]], na.rm = TRUE)
    
    scale_color_pca <- function() {
      ggplot2::scale_color_viridis_c(
        option = "viridis",
        limits = order_range,
        name = col
      )
    }
  } else {
    lvls <- sort(unique(c(before_df[[col]], after_df[[col]])))
    if (length(lvls) > length(cbPalette)) {
      stop("Too many groups for palette.")
    }
    
    cols <- stats::setNames(cbPalette[seq_along(lvls)], lvls)
    
    combined[[col]] <- factor(combined[[col]], levels = lvls)
    before_df[[col]] <- factor(before_df[[col]], levels = lvls)
    after_df[[col]] <- factor(after_df[[col]], levels = lvls)
    
    scale_color_pca <- function() {
      ggplot2::scale_color_manual(
        values = cols,
        name = col,
        drop = FALSE,
        na.translate = FALSE
      )
    }
  }
  
  big_font_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  p1 <- ggplot2::ggplot(before_df, ggplot2::aes(PC1, PC2, color = .data[[col]])) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = "Before",
      x = sprintf("PC1 (%.1f%%)", var_raw[1]),
      y = sprintf("PC2 (%.1f%%)", var_raw[2])
    ) +
    ggplot2::xlim(x_limits) +
    ggplot2::ylim(y_limits) +
    scale_color_pca() +
    big_font_theme +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 5, 10, 5)
    )
  
  p2 <- ggplot2::ggplot(after_df, ggplot2::aes(PC1, PC2, color = .data[[col]])) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = "After",
      x = sprintf("PC1 (%.1f%%)", var_cor[1]),
      y = sprintf("PC2 (%.1f%%)", var_cor[2])
    ) +
    ggplot2::xlim(x_limits) +
    ggplot2::ylim(y_limits) +
    scale_color_pca() +
    big_font_theme +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 5, 10, 5)
    )
  
  p_leg <- ggplot2::ggplot(before_df, ggplot2::aes(PC1, PC2, color = .data[[col]])) +
    ggplot2::geom_point(size = 2) +
    scale_color_pca() +
    ggplot2::guides(
      fill = "none",
      size = "none",
      shape = "none",
      alpha = "none",
      linetype = "none"
    ) +
    big_font_theme +
    ggplot2::theme(
      legend.position = "right",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    )
  
  leg <- cowplot::get_legend(p_leg)
  
  comb <- cowplot::plot_grid(
    p1, p2, leg,
    nrow = 1,
    rel_widths = c(1, 1, 0.3),
    labels = NULL,
    align = "hv",
    axis = "tblr"
  )
  
  cowplot::ggdraw() +
    cowplot::draw_label(
      paste("Comparison of PCA Before and After", compared_to),
      fontface = "bold",
      x = 0.5,
      y = 0.98,
      hjust = 0.5,
      vjust = 1,
      size = 14
    ) +
    cowplot::draw_plot(
      comb,
      x = 0,
      y = 0,
      width = 1,
      height = 0.93
    )
}


#' Build PCA loading plot from precomputed PCA results
#'
#' @param pca_pair list
#'   Output from \code{compute_pca_pair()}.
#' @param compared_to character
#'   Text appended to the title.
#' @param top_n integer
#'   Number of top loadings to display per PC.
#' @param label_width integer
#'   Width used for wrapped variable labels.
#'
#' @return ggplot
#'   Combined PCA loading plot.
#'
#' @noRd
plot_pca_loading_from_result <- function(
    pca_pair,
    compared_to,
    top_n = 5,
    label_width = 28
) {
  tidy_top <- function(loadings_df, label) {
    keep_cols <- intersect(c("PC1", "PC2"), names(loadings_df))
    if (length(keep_cols) == 0L) {
      stop("PCA loading data does not contain PC1 or PC2.")
    }
    
    df <- loadings_df[, c("variable", keep_cols), drop = FALSE]
    df <- tidyr::pivot_longer(
      df,
      cols = dplyr::all_of(keep_cols),
      names_to = "PC",
      values_to = "loading"
    )
    
    df <- dplyr::mutate(
      df,
      abs_loading = abs(.data$loading),
      sign = ifelse(.data$loading >= 0, "Positive", "Negative"),
      panel = label
    )
    
    df <- df |>
      dplyr::group_by(.data$PC) |>
      dplyr::slice_max(.data$abs_loading, n = top_n, with_ties = FALSE) |>
      dplyr::ungroup()
    
    df <- df |>
      dplyr::group_by(.data$panel, .data$PC) |>
      dplyr::arrange(dplyr::desc(.data$abs_loading), .by_group = TRUE) |>
      dplyr::mutate(
        variable_wrapped = stringr::str_wrap(.data$variable, width = label_width),
        variable_wrapped = factor(.data$variable_wrapped, levels = rev(unique(.data$variable_wrapped)))
      ) |>
      dplyr::ungroup()
    
    df
  }
  
  before_top <- tidy_top(pca_pair$before$loadings, "Before")
  after_top <- tidy_top(pca_pair$after$loadings, "After")
  
  lim <- max(c(before_top$abs_loading, after_top$abs_loading), na.rm = TRUE)
  y_limits <- c(-lim, lim)
  
  big_font_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  mk_plot <- function(df, title) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$variable_wrapped, y = .data$loading, fill = .data$sign)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::facet_wrap(~ PC, nrow = 2, scales = "free_y") +
      ggplot2::scale_y_continuous(
        limits = y_limits,
        expand = ggplot2::expansion(mult = c(0.02, 0.05))
      ) +
      ggplot2::scale_fill_manual(
        values = c("Positive" = "#4C9F50", "Negative" = "#BE0032"),
        drop = FALSE
      ) +
      ggplot2::labs(title = title, x = NULL, y = "Loading") +
      big_font_theme +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "none",
        plot.margin = ggplot2::margin(10, 5, 10, 5)
      )
  }
  
  p_before <- mk_plot(before_top, "Before")
  p_after <- mk_plot(after_top, "After")
  
  p_leg <- ggplot2::ggplot(before_top, ggplot2::aes(x = variable_wrapped, y = loading, fill = sign)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("Positive" = "#4C9F50", "Negative" = "#BE0032"),
      name = "Sign",
      drop = FALSE
    ) +
    ggplot2::guides(
      y = "none",
      x = "none",
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    big_font_theme +
    ggplot2::theme(legend.position = "right")
  
  leg <- cowplot::get_legend(p_leg)
  
  comb <- cowplot::plot_grid(
    p_before, p_after, leg,
    nrow = 1,
    rel_widths = c(0.7, 0.7, 0.22),
    labels = NULL,
    align = "hv",
    axis = "tblr"
  )
  
  cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("Top ", top_n, " Loadings for PC1 and PC2, Before vs After ", compared_to),
      fontface = "bold",
      x = 0.5,
      y = 0.98,
      hjust = 0.5,
      vjust = 1,
      size = 14
    ) +
    cowplot::draw_plot(comb, x = 0, y = 0, width = 1, height = 0.93)
}


#' Backward-compatible PCA score plotting wrapper
#'
#' @param p list
#'   Parameter list.
#' @param before data.frame
#'   Before data set.
#' @param after data.frame
#'   After data set.
#' @param compared_to character
#'   Text appended to the title.
#'
#' @return ggplot
#'
#' @noRd
plot_pca <- function(p, before, after, compared_to) {
  pca_pair <- compute_pca_pair(
    before = before,
    after = after,
    p = p,
    before_label = "Before",
    after_label = "After"
  )
  
  plot_pca_from_result(
    p = p,
    pca_pair = pca_pair,
    compared_to = compared_to
  )
}


#' Backward-compatible PCA loading plotting wrapper
#'
#' @param p list
#'   Parameter list.
#' @param before data.frame
#'   Before data set.
#' @param after data.frame
#'   After data set.
#' @param compared_to character
#'   Text appended to the title.
#' @param top_n integer
#'   Number of top loadings to display per PC.
#' @param label_width integer
#'   Width used for wrapped variable labels.
#'
#' @return ggplot
#'
#' @noRd
plot_pca_loading <- function(
    p,
    before,
    after,
    compared_to,
    top_n = 5,
    label_width = 28
) {
  pca_pair <- compute_pca_pair(
    before = before,
    after = after,
    p = p,
    before_label = "Before",
    after_label = "After"
  )
  
  plot_pca_loading_from_result(
    pca_pair = pca_pair,
    compared_to = compared_to,
    top_n = top_n,
    label_width = label_width
  )
}


#' Add summary columns to a loading table from precomputed PCA result
#'
#' @param pca_res list
#'   Output from \code{compute_single_pca()}.
#' @param dataset_label character
#'   Label identifying the dataset.
#'
#' @return list
#'   A list containing loadings and explained variance data frames.
#'
#' @noRd
compute_pca_loadings_table_from_result <- function(pca_res, dataset_label) {
  loadings_df <- pca_res$loadings
  
  pc_cols <- grep("^PC\\d+$", names(loadings_df), value = TRUE)
  
  if ("PC1" %in% pc_cols) {
    loadings_df$abs_PC1 <- abs(loadings_df$PC1)
  }
  if ("PC2" %in% pc_cols) {
    loadings_df$abs_PC2 <- abs(loadings_df$PC2)
  }
  
  if (length(pc_cols) > 0L) {
    loadings_df$max_abs_loading <- apply(
      abs(loadings_df[, pc_cols, drop = FALSE]),
      MARGIN = 1,
      FUN = max,
      na.rm = TRUE
    )
  } else {
    loadings_df$max_abs_loading <- NA_real_
  }
  
  loadings_df$dataset <- dataset_label
  loadings_df <- loadings_df[, c("dataset", setdiff(names(loadings_df), "dataset")), drop = FALSE]
  
  explained_variance_df <- pca_res$explained_variance
  explained_variance_df$dataset <- dataset_label
  explained_variance_df <- explained_variance_df[, c("dataset", setdiff(names(explained_variance_df), "dataset")), drop = FALSE]
  
  list(
    loadings = loadings_df,
    explained_variance = explained_variance_df
  )
}


#' Compute PCA loading/export tables from a raw data frame
#'
#' @param df data.frame
#'   Input data frame containing metadata and metabolite columns.
#' @param p list
#'   Parameter list containing imputation settings.
#' @param dataset_label character
#'   Label identifying the dataset.
#' @param meta_cols character
#'   Metadata columns to exclude.
#'
#' @return list
#'   A list containing loadings and explained variance data frames.
#'
#' @noRd
compute_pca_loadings_table <- function(
    df,
    p,
    dataset_label,
    meta_cols = c("sample", "batch", "class", "order")
) {
  metab_cols <- setdiff(names(df), meta_cols)
  
  pca_res <- compute_single_pca(
    df = df,
    p = p,
    metab_cols = metab_cols,
    meta_cols = meta_cols
  )
  
  compute_pca_loadings_table_from_result(
    pca_res = pca_res,
    dataset_label = dataset_label
  )
}


#' Get data frames used for PCA export
#'
#' @param p list
#'   Parameter list.
#' @param d list
#'   Data object containing filtered/corrected/transformed data.
#'
#' @return named list
#'   Named list of data frames to export.
#'
#' @noRd
get_pca_export_datasets <- function(p, d) {
  datasets <- list()
  
  datasets$raw_data <- d$filtered$df
  
  datasets$corrected_data <- if (isTRUE(p$remove_imputed)) {
    d$filtered_corrected$df_mv
  } else {
    d$filtered_corrected$df_no_mv
  }
  
  if (!identical(p$transform, "none")) {
    datasets$transformed_data <- if (isTRUE(p$remove_imputed)) {
      d$transformed$df_mv
    } else {
      d$transformed$df_no_mv
    }
  }
  
  datasets
}


#' Compute PCA results for all exportable PCA datasets
#'
#' @param p list
#'   Parameter list.
#' @param d list
#'   Data object containing filtered/corrected/transformed data.
#' @param meta_cols character
#'   Metadata columns to retain.
#'
#' @return named list
#'   Named list of PCA results.
#'
#' @noRd
compute_all_pca_export_results <- function(
    p,
    d,
    meta_cols = c("sample", "batch", "class", "order")
) {
  datasets <- get_pca_export_datasets(p, d)
  
  lapply(datasets, function(df) {
    metab_cols <- setdiff(names(df), meta_cols)
    compute_single_pca(
      df = df,
      p = p,
      metab_cols = metab_cols,
      meta_cols = meta_cols
    )
  })
}


#' Export PCA loadings to Excel
#'
#' Creates an Excel workbook containing PCA loadings for raw, corrected, and
#' transformed data. Each dataset gets its own worksheet. An additional
#' explained-variance worksheet is also included.
#'
#' @param p list
#'   Parameter list.
#' @param d list
#'   Data object containing analysis data frames.
#' @param pca_dir character
#'   Directory where the workbook should be saved.
#' @param file_name character
#'   Output workbook file name.
#' @param pca_results named list or NULL
#'   Optional precomputed PCA results from
#'   \code{compute_all_pca_export_results()}.
#'
#' @return character
#'   Normalized path to the written workbook.
#'
#' @noRd
export_pca_loadings_xlsx <- function(
    p,
    d,
    pca_dir,
    file_name = "pca_loadings.xlsx",
    pca_results = NULL
) {
  .require_pkg("openxlsx", "write PCA loadings workbook")
  
  if (is.null(pca_results)) {
    pca_results <- compute_all_pca_export_results(p = p, d = d)
  }
  
  wb <- openxlsx::createWorkbook()
  
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "bottom"
  )
  
  num_style <- openxlsx::createStyle(numFmt = "0.0000")
  
  ev_all <- list()
  
  for (nm in names(pca_results)) {
    res <- compute_pca_loadings_table_from_result(
      pca_res = pca_results[[nm]],
      dataset_label = nm
    )
    
    sheet_name <- switch(
      nm,
      raw_data = "Raw PCA Loadings",
      corrected_data = "Corrected PCA Loadings",
      transformed_data = "Transformed PCA Loadings",
      nm
    )
    
    sheet_name <- substr(sheet_name, 1L, 31L)
    
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeDataTable(
      wb,
      sheet = sheet_name,
      x = res$loadings,
      withFilter = TRUE,
      tableStyle = "TableStyleMedium2"
    )
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = header_style,
      rows = 1,
      cols = seq_len(ncol(res$loadings)),
      gridExpand = TRUE,
      stack = TRUE
    )
    
    numeric_cols <- which(vapply(res$loadings, is.numeric, logical(1)))
    if (length(numeric_cols) > 0L) {
      openxlsx::addStyle(
        wb,
        sheet = sheet_name,
        style = num_style,
        rows = 2:(nrow(res$loadings) + 1L),
        cols = numeric_cols,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
    
    openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
    openxlsx::setColWidths(
      wb,
      sheet = sheet_name,
      cols = 1:ncol(res$loadings),
      widths = "auto"
    )
    
    ev_all[[nm]] <- res$explained_variance
  }
  
  explained_variance_df <- do.call(rbind, ev_all)
  
  openxlsx::addWorksheet(wb, "Explained Variance")
  openxlsx::writeDataTable(
    wb,
    sheet = "Explained Variance",
    x = explained_variance_df,
    withFilter = TRUE,
    tableStyle = "TableStyleMedium2"
  )
  
  openxlsx::addStyle(
    wb,
    sheet = "Explained Variance",
    style = header_style,
    rows = 1,
    cols = seq_len(ncol(explained_variance_df)),
    gridExpand = TRUE,
    stack = TRUE
  )
  
  numeric_cols_ev <- which(vapply(explained_variance_df, is.numeric, logical(1)))
  if (length(numeric_cols_ev) > 0L) {
    openxlsx::addStyle(
      wb,
      sheet = "Explained Variance",
      style = num_style,
      rows = 2:(nrow(explained_variance_df) + 1L),
      cols = numeric_cols_ev,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  openxlsx::freezePane(wb, sheet = "Explained Variance", firstRow = TRUE)
  openxlsx::setColWidths(
    wb,
    sheet = "Explained Variance",
    cols = 1:ncol(explained_variance_df),
    widths = "auto"
  )
  
  out_path <- file.path(pca_dir, file_name)
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  
  normalizePath(out_path, winslash = "/", mustWork = TRUE)
}


#' Get before/after data for a PCA comparison option
#'
#' @param p list
#'   Parameter list.
#' @param d list
#'   Data object.
#' @param pca_compare character
#'   PCA comparison option.
#'
#' @return list
#'   A list with before, after, and compared_to.
#'
#' @noRd
get_pca_compare_data <- function(p, d, pca_compare) {
  before <- d$filtered$df
  
  after <- switch(
    pca_compare,
    filtered_cor_data = if (isTRUE(p$remove_imputed)) d$filtered_corrected$df_mv else d$filtered_corrected$df_no_mv,
    transformed_cor_data = if (isTRUE(p$remove_imputed)) d$transformed$df_mv else d$transformed$df_no_mv,
    stop(sprintf("Unsupported pca_compare value: %s", pca_compare))
  )
  
  compared_to <- switch(
    pca_compare,
    filtered_cor_data = "Corrected Data",
    transformed_cor_data = "Transformed Data",
    pca_compare
  )
  
  list(
    before = before,
    after = after,
    compared_to = compared_to
  )
}


#' Make all PCA score and loading plots
#'
#' Computes each PCA comparison once and reuses the result for score plots
#' and loading plots.
#'
#' @param p list
#'   Parameter list.
#' @param d list
#'   Data object.
#'
#' @return list
#'   A list containing PCA plots, names, loading plots, loading names, and
#'   paired PCA results.
#'
#' @noRd
make_all_pca_plots <- function(p, d) {
  build_name <- function(compare, color) {
    sprintf("pca_%s_%s", compare, color)
  }
  
  specs <- expand.grid(
    color_col = c("batch", "class", "order"),
    pca_compare = "filtered_cor_data",
    stringsAsFactors = FALSE
  )
  
  if (!identical(p$transform, "none")) {
    specs_trans <- specs
    specs_trans$pca_compare <- "transformed_cor_data"
    specs <- rbind(specs, specs_trans)
  }
  
  unique_compares <- unique(specs$pca_compare)
  
  pca_pairs <- stats::setNames(vector("list", length(unique_compares)), unique_compares)
  for (cmp in unique_compares) {
    compare_data <- get_pca_compare_data(p = p, d = d, pca_compare = cmp)
    pca_pairs[[cmp]] <- compute_pca_pair(
      before = compare_data$before,
      after = compare_data$after,
      p = p,
      before_label = "Before",
      after_label = "After"
    )
    pca_pairs[[cmp]]$compared_to <- compare_data$compared_to
  }
  
  pca_plots <- vector("list", nrow(specs))
  plot_names <- character(nrow(specs))
  
  for (i in seq_len(nrow(specs))) {
    temp_params <- p
    temp_params$color_col <- specs$color_col[i]
    temp_params$pca_compare <- specs$pca_compare[i]
    
    pca_plots[[i]] <- plot_pca_from_result(
      p = temp_params,
      pca_pair = pca_pairs[[temp_params$pca_compare]],
      compared_to = pca_pairs[[temp_params$pca_compare]]$compared_to
    )
    
    plot_names[i] <- build_name(temp_params$pca_compare, temp_params$color_col)
  }
  
  loading_compares <- unique(specs$pca_compare)
  pca_loading_plots <- vector("list", length(loading_compares))
  loading_plot_names <- character(length(loading_compares))
  
  for (i in seq_along(loading_compares)) {
    cmp <- loading_compares[i]
    pca_loading_plots[[i]] <- plot_pca_loading_from_result(
      pca_pair = pca_pairs[[cmp]],
      compared_to = pca_pairs[[cmp]]$compared_to
    )
    loading_plot_names[i] <- sprintf("pca_loadings_%s", cmp)
  }
  
  list(
    pca_plots = pca_plots,
    plot_names = plot_names,
    pca_loading_plots = pca_loading_plots,
    loading_plot_names = loading_plot_names,
    pca_pairs = pca_pairs
  )
}