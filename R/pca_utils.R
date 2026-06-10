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
#' Validate PCA plotting metadata
#'
#' Ensures the metadata contains a sample column, has unique samples, and
#' includes all samples present in the PCA data frame.
#'
#' @param df data.frame
#'   PCA input data frame.
#' @param meta_df data.frame
#'   Metadata data frame used for plotting.
#' @param sample_col character
#'   Sample identifier column.
#'
#' @return data.frame
#'   A validated metadata data frame.
#'
#' @noRd
validate_pca_meta_df <- function(df, meta_df, sample_col = "sample") {
  if (is.null(meta_df)) {
    stop("`meta_df` cannot be NULL.")
  }

  if (!sample_col %in% names(df)) {
    stop(sprintf("PCA data frame must contain '%s'.", sample_col))
  }

  if (!sample_col %in% names(meta_df)) {
    stop(sprintf("`meta_df` must contain '%s'.", sample_col))
  }

  if (anyDuplicated(meta_df[[sample_col]]) > 0L) {
    dupes <- unique(meta_df[[sample_col]][duplicated(meta_df[[sample_col]])])
    stop(
      sprintf(
        "`meta_df` contains duplicate sample identifiers. Examples: %s",
        paste(utils::head(dupes, 10L), collapse = ", ")
      )
    )
  }

  missing_meta <- setdiff(df[[sample_col]], meta_df[[sample_col]])
  if (length(missing_meta) > 0L) {
    stop(
      sprintf(
        "Some PCA samples are missing from `meta_df`. Examples: %s",
        paste(utils::head(missing_meta, 10L), collapse = ", ")
      )
    )
  }

  meta_df
}


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
#' validates that all PCA inputs are numeric, imputes missing values with KNN
#' if needed, and verifies that no missing values remain before PCA.
#'
#' PCA uses KNN imputation regardless of the user-selected QC/sample imputation
#' settings because `stats::prcomp()` cannot handle missing values.
#'
#' @param df data.frame
#'   Input data frame containing metadata and metabolite columns.
#' @param p list
#'   Parameter list. Included for compatibility with existing PCA calls.
#'   User-selected imputation methods are intentionally ignored for PCA.
#' @param metab_cols character
#'   Metabolite columns to include in PCA.
#'
#' @return data.frame
#'   Numeric data frame containing only metabolite columns, with missing values
#'   imputed by KNN when needed.
#'
#' @noRd
prep_pca_matrix <- function(df, p, metab_cols) {
  if (length(metab_cols) == 0L) {
    stop("No metabolite columns available for PCA.")
  }

  missing_metab_cols <- setdiff(metab_cols, names(df))
  if (length(missing_metab_cols) > 0L) {
    stop(
      sprintf(
        "These PCA metabolite columns are missing from `df`: %s",
        paste(utils::head(missing_metab_cols, 20L), collapse = ", ")
      )
    )
  }

  out <- df[, metab_cols, drop = FALSE]

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

  all_missing_cols <- names(out)[
    vapply(out, function(x) all(is.na(x)), logical(1))
  ]

  if (length(all_missing_cols) > 0L) {
    stop(
      sprintf(
        paste(
          "KNN imputation cannot be used for PCA because these metabolite",
          "columns are entirely missing: %s"
        ),
        paste(utils::head(all_missing_cols, 20L), collapse = ", ")
      )
    )
  }

  if (anyNA(out)) {
    results <- impute_missing(
      df = df,
      metab_cols = metab_cols,
      qcImputeM = "KNN",
      samImputeM = "KNN"
    )

    out <- results$df[, metab_cols, drop = FALSE]

    if (anyNA(out)) {
      still_missing <- names(out)[
        vapply(out, function(x) any(is.na(x)), logical(1))
      ]

      stop(
        sprintf(
          paste(
            "KNN imputation was applied before PCA, but missing values remain",
            "in these metabolite columns: %s"
          ),
          paste(utils::head(still_missing, 20L), collapse = ", ")
        )
      )
    }
  }

  out
}


#' Compute PCA for a single data frame
#'
#' Fits PCA once and returns the model, sample scores, variable loadings,
#' and explained variance.
#'
#' @param df data.frame
#'   Input data frame containing metabolite columns and a sample column.
#' @param p list
#'   Parameter list containing imputation settings.
#' @param metab_cols character
#'   Metabolite columns to include in PCA.
#' @param meta_cols character
#'   Metadata columns to retain in the score output.
#' @param meta_df data.frame or NULL
#'   Optional external metadata data frame used for plotting. Must contain
#'   at least the sample column and any requested metadata columns.
#' @param sample_col character
#'   Sample identifier column used to align metadata.
#'
#' @return list
#'   A list containing PCA fit, scores, loadings, explained variance,
#'   and metabolite columns used.
#'
#' @noRd
compute_single_pca <- function(
  df,
  p,
  metab_cols,
  meta_cols = c("sample", "batch", "class", "order"),
  meta_df = NULL,
  sample_col = "sample"
) {
  x <- prep_pca_matrix(df = df, p = p, metab_cols = metab_cols)
  fit <- stats::prcomp(x, center = TRUE, scale. = TRUE)

  if (is.null(meta_df)) {
    meta_source <- df
  } else {
    if (!sample_col %in% names(df)) {
      stop(sprintf("`df` must contain '%s'.", sample_col))
    }
    if (!sample_col %in% names(meta_df)) {
      stop(sprintf("`meta_df` must contain '%s'.", sample_col))
    }
    if (anyDuplicated(meta_df[[sample_col]]) > 0L) {
      stop("`meta_df` contains duplicate sample values.")
    }

    idx <- match(df[[sample_col]], meta_df[[sample_col]])
    if (anyNA(idx)) {
      missing_samples <- df[[sample_col]][is.na(idx)]
      stop(
        sprintf(
          "Some PCA samples are missing from `meta_df`. Examples: %s",
          paste(utils::head(missing_samples, 10L), collapse = ", ")
        )
      )
    }

    meta_source <- meta_df[idx, , drop = FALSE]
  }

  available_meta <- intersect(meta_cols, names(meta_source))

  if (!sample_col %in% available_meta && sample_col %in% names(meta_source)) {
    available_meta <- c(sample_col, available_meta)
  }

  scores_df <- as.data.frame(fit$x)
  scores_df <- dplyr::bind_cols(scores_df, meta_source[, available_meta, drop = FALSE])

  loadings_df <- as.data.frame(fit$rotation)
  loadings_df$variable <- rownames(loadings_df)
  rownames(loadings_df) <- NULL
  loadings_df <- loadings_df[, c("variable", setdiff(names(loadings_df), "variable")), drop = FALSE]

  explained_var <- (fit$sdev^2) / sum(fit$sdev^2)
  explained_variance_df <- data.frame(
    PC = paste0("PC", seq_along(explained_var)),
    explained_variance = as.numeric(explained_var),
    cumulative_explained_variance = cumsum(as.numeric(explained_var))
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
#' @param meta_df data.frame or NULL
#'   Optional external metadata data frame used for coloring and labeling.
#' @param sample_col character
#'   Sample identifier column.
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
  meta_cols = c("sample", "batch", "class", "order"),
  meta_df = NULL,
  sample_col = "sample"
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
    meta_cols = meta_cols,
    meta_df = meta_df,
    sample_col = sample_col
  )

  after_res <- compute_single_pca(
    df = after,
    p = p,
    metab_cols = metab_cols,
    meta_cols = meta_cols,
    meta_df = meta_df,
    sample_col = sample_col
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
#'   Parameter list. Uses `p$color_col` for point colors and optionally
#'   `p$shape_col` for point shapes.
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

  pc1_range <- range(combined$PC1, na.rm = TRUE)
  pc2_range <- range(combined$PC2, na.rm = TRUE)
  max_abs <- max(abs(c(pc1_range, pc2_range)), na.rm = TRUE)
  axis_limits <- c(-max_abs, max_abs)

  var_raw <- 100 * pca_pair$before$explained_variance$explained_variance[1:2]
  var_cor <- 100 * pca_pair$after$explained_variance$explained_variance[1:2]

  cbPalette <- c(
    "#F3C300", "#875692", "#ee7733", "#A1CAF1", "#BE0032",
    "#C2B280", "#555555", "#008856", "#E68FAC", "#0067A5",
    "#F99379", "#332288", "#F6A600", "#B3446C", "#DCD300",
    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
    "#bbbbbb", "#000000", "#33bbee", "#ccddaa", "#225555"
  )

  pca_shapes <- c(
    16, 17, 15, 18, 1, 2, 0, 5, 6, 3, 7, 8,
    9, 10, 11, 12, 13, 14, 4, 19, 20, 21, 22, 23, 24, 25
  )

  col <- as.character(p$color_col %||% "class")
  shape_col <- as.character(p$shape_col %||% "none")

  use_shape <- !is.null(shape_col) &&
    length(shape_col) == 1L &&
    !is.na(shape_col) &&
    nzchar(shape_col) &&
    !identical(shape_col, "none")

  if (!col %in% names(before_df) || !col %in% names(after_df)) {
    stop(sprintf("Column '%s' not found in PCA score data.", col))
  }

  if (isTRUE(use_shape)) {
    if (!shape_col %in% names(before_df) || !shape_col %in% names(after_df)) {
      stop(sprintf("Shape column '%s' not found in PCA score data.", shape_col))
    }
  }

  is_numeric_like <- function(x) {
    is.numeric(x) || is.integer(x)
  }

  before_is_numeric <- is_numeric_like(before_df[[col]])
  after_is_numeric <- is_numeric_like(after_df[[col]])

  if (before_is_numeric != after_is_numeric) {
    stop(sprintf(
      "Column '%s' is not of the same type in before/after PCA score data.",
      col
    ))
  }

  use_gradient <- before_is_numeric && after_is_numeric
  legend_ncol <- 1L
  legend_rel_width <- 0.32

  if (use_gradient) {
    combined[[col]] <- as.numeric(combined[[col]])
    before_df[[col]] <- as.numeric(before_df[[col]])
    after_df[[col]] <- as.numeric(after_df[[col]])

    color_range <- range(combined[[col]], na.rm = TRUE)
    legend_rel_width <- 0.32

    scale_color_pca <- function() {
      ggplot2::scale_color_viridis_c(
        option = "viridis",
        limits = color_range,
        name = col,
        na.value = "grey80"
      )
    }
  } else {
    lvls <- sort(unique(as.character(c(before_df[[col]], after_df[[col]]))))

    if (length(lvls) > length(cbPalette)) {
      stop("Too many groups for palette.")
    }

    cols <- stats::setNames(cbPalette[seq_along(lvls)], lvls)

    before_df[[col]] <- factor(as.character(before_df[[col]]), levels = lvls)
    after_df[[col]] <- factor(as.character(after_df[[col]]), levels = lvls)

    if (length(lvls) > 12L) {
      legend_ncol <- 2L
      legend_rel_width <- 0.65
    }

    scale_color_pca <- function() {
      ggplot2::scale_color_manual(
        values = cols,
        name = col,
        drop = FALSE,
        na.translate = FALSE
      )
    }
  }

  shape_aes_col <- ".pca_shape_group"
  shape_legend_ncol <- 1L

  if (isTRUE(use_shape)) {
    before_shape_is_numeric <- is_numeric_like(before_df[[shape_col]])
    after_shape_is_numeric <- is_numeric_like(after_df[[shape_col]])

    if (before_shape_is_numeric != after_shape_is_numeric) {
      stop(sprintf(
        "Shape column '%s' is not of the same type in before/after PCA score data.",
        shape_col
      ))
    }

    if (before_shape_is_numeric && after_shape_is_numeric) {
      shape_lvls <- sort(unique(as.numeric(c(
        before_df[[shape_col]],
        after_df[[shape_col]]
      ))))
      shape_lvls <- as.character(shape_lvls)
    } else {
      shape_lvls <- sort(unique(as.character(c(
        before_df[[shape_col]],
        after_df[[shape_col]]
      ))))
    }

    shape_lvls <- shape_lvls[!is.na(shape_lvls)]

    if (length(shape_lvls) == 0L) {
      stop(sprintf("Shape column '%s' contains no non-missing values.", shape_col))
    }

    if (length(shape_lvls) > length(pca_shapes)) {
      stop(sprintf(
        paste(
          "Too many groups for point shapes in column '%s'.",
          "Found %d groups, but only %d shapes are available."
        ),
        shape_col,
        length(shape_lvls),
        length(pca_shapes)
      ))
    }

    shape_values <- stats::setNames(pca_shapes[seq_along(shape_lvls)], shape_lvls)

    before_df[[shape_aes_col]] <- factor(
      as.character(before_df[[shape_col]]),
      levels = shape_lvls
    )

    after_df[[shape_aes_col]] <- factor(
      as.character(after_df[[shape_col]]),
      levels = shape_lvls
    )

    if (length(shape_lvls) > 8L) {
      shape_legend_ncol <- 2L
      legend_rel_width <- max(legend_rel_width, 0.55)
    } else {
      legend_rel_width <- max(legend_rel_width, 0.38)
    }

    scale_shape_pca <- function() {
      ggplot2::scale_shape_manual(
        values = shape_values,
        name = shape_col,
        drop = FALSE,
        na.translate = FALSE
      )
    }
  } else {
    scale_shape_pca <- function() {
      NULL
    }
  }

  legend_guides <- function() {
    if (isTRUE(use_shape)) {
      ggplot2::guides(
        color = if (use_gradient) {
          ggplot2::guide_colorbar(title.position = "top")
        } else {
          ggplot2::guide_legend(
            ncol = legend_ncol,
            byrow = TRUE,
            title.position = "top",
            override.aes = list(size = 3)
          )
        },
        shape = ggplot2::guide_legend(
          ncol = shape_legend_ncol,
          byrow = TRUE,
          title.position = "top",
          override.aes = list(size = 3)
        ),
        fill = "none",
        size = "none",
        alpha = "none",
        linetype = "none"
      )
    } else {
      ggplot2::guides(
        color = if (use_gradient) {
          ggplot2::guide_colorbar(title.position = "top")
        } else {
          ggplot2::guide_legend(
            ncol = legend_ncol,
            byrow = TRUE,
            title.position = "top"
          )
        },
        fill = "none",
        size = "none",
        shape = "none",
        alpha = "none",
        linetype = "none"
      )
    }
  }

  point_mapping <- if (isTRUE(use_shape)) {
    ggplot2::aes(
      x = .data$PC1,
      y = .data$PC2,
      color = .data[[col]],
      shape = .data[[shape_aes_col]]
    )
  } else {
    ggplot2::aes(
      x = .data$PC1,
      y = .data$PC2,
      color = .data[[col]]
    )
  }

  point_size <- if (isTRUE(use_shape)) 2.5 else 2

  big_font_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )

  panel_theme <- ggplot2::theme(
    legend.position = "none",
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = ggplot2::margin(10, 5, 10, 5)
  )

  p1 <- ggplot2::ggplot(before_df, point_mapping) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::labs(
      title = "Before",
      x = sprintf("PC1 (%.1f%%)", var_raw[1]),
      y = sprintf("PC2 (%.1f%%)", var_raw[2])
    ) +
    scale_color_pca() +
    scale_shape_pca() +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = axis_limits,
      ylim = axis_limits,
      expand = TRUE,
      clip = "on"
    ) +
    big_font_theme +
    panel_theme

  p2 <- ggplot2::ggplot(after_df, point_mapping) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::labs(
      title = "After",
      x = sprintf("PC1 (%.1f%%)", var_cor[1]),
      y = sprintf("PC2 (%.1f%%)", var_cor[2])
    ) +
    scale_color_pca() +
    scale_shape_pca() +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = axis_limits,
      ylim = axis_limits,
      expand = TRUE,
      clip = "on"
    ) +
    big_font_theme +
    panel_theme

  p_leg <- ggplot2::ggplot(before_df, point_mapping) +
    ggplot2::geom_point(size = point_size) +
    scale_color_pca() +
    scale_shape_pca() +
    legend_guides() +
    big_font_theme +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.key.width = grid::unit(0.45, "cm")
    )

  leg <- cowplot::get_legend(p_leg)

  comb <- cowplot::plot_grid(
    p1,
    p2,
    leg,
    nrow = 1,
    rel_widths = c(1, 1, legend_rel_width),
    align = "h",
    axis = "tb"
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
      ggplot2::facet_wrap(~PC, nrow = 2, scales = "free_y") +
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
  bold <- openxlsx::createStyle(textDecoration = "Bold")
  note <- openxlsx::createStyle(
    wrapText = TRUE,
    valign = "top",
    fgFill = "#f8cbad"
  )

  num_style <- openxlsx::createStyle(numFmt = "0.0000")

  ev_all <- list()

  # leave space for title + description + blank row
  data_start_row <- 5L

  for (nm in names(pca_results)) {
    res <- compute_pca_loadings_table_from_result(
      pca_res = pca_results[[nm]],
      dataset_label = nm
    )

    sheet_name <- switch(nm,
      raw_data = "Raw PCA Loadings",
      corrected_data = "Corrected PCA Loadings",
      transformed_data = "Transformed PCA Loadings",
      nm
    )

    sheet_name <- substr(sheet_name, 1L, 31L)

    description_text <- paste(
      "PCA loadings describe how strongly each metabolite contributes to each principal component (PC).",
      "Larger absolute loading values indicate stronger influence on that PC.",
      "The sign indicates direction along the component, but the magnitude is usually the more important quantity when identifying influential metabolites.",
      "Columns abs_PC1 and abs_PC2 give the absolute loading magnitude for PC1 and PC2, and max_abs_loading gives the largest absolute loading across all PCs for each metabolite."
    )

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = res$loadings,
      startRow = data_start_row,
      startCol = 1,
      rowNames = FALSE,
      headerStyle = bold,
      withFilter = TRUE
    )

    numeric_cols <- which(vapply(res$loadings, is.numeric, logical(1)))
    if (length(numeric_cols) > 0L) {
      openxlsx::addStyle(
        wb,
        sheet = sheet_name,
        style = num_style,
        rows = (data_start_row + 1L):(data_start_row + nrow(res$loadings)),
        cols = numeric_cols,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
    openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1, widths = 25)
    openxlsx::setColWidths(wb, sheet = sheet_name, cols = 2:ncol(res$loadings), widths = "auto")

    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = description_text,
      startRow = 1,
      startCol = 1,
      colNames = FALSE,
      rowNames = FALSE
    )
    openxlsx::mergeCells(wb, sheet = sheet_name, cols = 1:8, rows = 1)
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, sheet = sheet_name, rows = 1, heights = 60)

    ev_all[[nm]] <- res$explained_variance
  }

  explained_variance_df <- do.call(rbind, ev_all)

  ev_sheet_name <- "Explained Variance"
  ev_description_text <- paste(
    "Explained variance gives the proportion of total variance captured by each principal component (PC).",
    "Higher explained_variance means that PC summarizes more of the structure in the data.",
    "cumulative_explained_variance shows the running total across PCs and is useful for assessing how many components are needed to represent the dataset."
  )

  openxlsx::addWorksheet(wb, ev_sheet_name)
  openxlsx::writeData(
    wb,
    sheet = ev_sheet_name,
    x = explained_variance_df,
    startRow = data_start_row,
    startCol = 1,
    rowNames = FALSE,
    headerStyle = bold,
    withFilter = TRUE
  )

  numeric_cols_ev <- which(vapply(explained_variance_df, is.numeric, logical(1)))
  if (length(numeric_cols_ev) > 0L) {
    openxlsx::addStyle(
      wb,
      sheet = ev_sheet_name,
      style = num_style,
      rows = (data_start_row + 1L):(data_start_row + nrow(explained_variance_df)),
      cols = numeric_cols_ev,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  openxlsx::setColWidths(wb, sheet = ev_sheet_name, cols = 1, widths = 25)
  openxlsx::setColWidths(
    wb,
    sheet = ev_sheet_name,
    cols = 2:ncol(explained_variance_df),
    widths = "auto"
  )
  openxlsx::writeData(
    wb,
    sheet = ev_sheet_name,
    x = ev_description_text,
    startRow = 1,
    startCol = 1,
    colNames = FALSE,
    rowNames = FALSE
  )

  openxlsx::mergeCells(wb, sheet = ev_sheet_name, cols = 1:8, rows = 1)
  openxlsx::addStyle(
    wb,
    sheet = ev_sheet_name,
    style = note,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  openxlsx::setRowHeights(wb, sheet = ev_sheet_name, rows = 1, heights = 60)

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

  after <- switch(pca_compare,
    filtered_cor_data = if (isTRUE(p$remove_imputed)) d$filtered_corrected$df_mv else d$filtered_corrected$df_no_mv,
    transformed_cor_data = if (isTRUE(p$remove_imputed)) d$transformed$df_mv else d$transformed$df_no_mv,
    stop(sprintf("Unsupported pca_compare value: %s", pca_compare))
  )

  compared_to <- switch(pca_compare,
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
#' @param p list
#'   Parameter list.
#' @param d list
#'   Data object.
#' @param meta_df data.frame or NULL
#'   Optional external metadata data frame for coloring PCA scores.
#'
#' @return list
#'   A list containing PCA plots, names, loading plots, loading names, and
#'   paired PCA results.
#'
#' @noRd
make_all_pca_plots <- function(p, d, meta_df = NULL) {
  build_name <- function(compare, color, shape) {
    sprintf("pca_%s_%s_%s", compare, color, shape)
  }

  if (is.null(meta_df)) {
    color_choices <- c("batch", "class", "order")
    shape_choices <- c("batch", "class")
    meta_cols <- c("sample", "batch", "class", "order")
  } else {
    color_choices <- setdiff(names(meta_df), "sample")
    shape_choices <- setdiff(names(meta_df), c("sample", "order"))
    meta_cols <- c("sample", color_choices)
  }

  specs <- expand.grid(
    color_col = color_choices,
    shape_col = shape_choices,
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
      after_label = "After",
      meta_cols = meta_cols,
      meta_df = meta_df,
      sample_col = "sample"
    )
    pca_pairs[[cmp]]$compared_to <- compare_data$compared_to
  }

  pca_plots <- vector("list", nrow(specs))
  plot_names <- character(nrow(specs))

  for (i in seq_len(nrow(specs))) {
    temp_params <- p
    temp_params$color_col <- specs$color_col[i]
    temp_params$shape_col <- specs$shape_col[i]
    temp_params$pca_compare <- specs$pca_compare[i]

    pca_plots[[i]] <- plot_pca_from_result(
      p = temp_params,
      pca_pair = pca_pairs[[temp_params$pca_compare]],
      compared_to = pca_pairs[[temp_params$pca_compare]]$compared_to
    )

    plot_names[i] <- build_name(temp_params$pca_compare, temp_params$color_col, temp_params$shape_col)
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
