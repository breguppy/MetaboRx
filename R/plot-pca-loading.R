#' PCA loading plots (PC1 & PC2), Before vs After
#' @keywords internal
#' @noRd
plot_pca_loading <- function(p, before, after, compared_to,
                             top_n = 5, label_width = 28) {
  # Overlapping metabolite columns
  meta_cols   <- c("sample", "batch", "class", "order")
  metab_cols  <- intersect(setdiff(names(before), meta_cols),
                           setdiff(names(after),  meta_cols))
  
  # Impute 'after' like in plot_pca()
  if (any(is.na(after[ , metab_cols, drop = FALSE]))) {
    results <- impute_missing(after, metab_cols, p$qcImputeM, p$samImputeM)
    after_df_only <- results$df
  } else {
    after_df_only <- after
  }
  
  # PCA
  before_p <- stats::prcomp(before[, metab_cols], center = TRUE, scale. = TRUE)
  after_p  <- stats::prcomp(after_df_only[, metab_cols], center = TRUE, scale. = TRUE)
  
  # Helper to tidy top loadings
  tidy_top <- function(rot, label) {
    df <- tibble::as_tibble(rot[, 1:2, drop = FALSE], rownames = "variable")
    df <- tidyr::pivot_longer(df, cols = c("PC1", "PC2"),
                              names_to = "PC", values_to = "loading")
    df <- dplyr::mutate(df,
                        abs_loading = abs(.data$loading),
                        sign        = ifelse(.data$loading >= 0, "Positive", "Negative"),
                        panel       = label)
    # top_n per PC
    df <- df |>
      dplyr::group_by(.data$PC) |>
      dplyr::slice_max(.data$abs_loading, n = top_n, with_ties = FALSE) |>
      dplyr::ungroup()
    
    # wrap labels and keep order by magnitude within each PC
    df <- df |>
      dplyr::group_by(.data$panel, .data$PC) |>
      dplyr::arrange(dplyr::desc(.data$abs_loading), .by_group = TRUE) |>
      dplyr::mutate(variable_wrapped = stringr::str_wrap(.data$variable, width = label_width),
                    variable_wrapped = factor(.data$variable_wrapped, levels = rev(unique(.data$variable_wrapped)))) |>
      dplyr::ungroup()
    df
  }
  
  before_top <- tidy_top(before_p$rotation, "Before")
  after_top  <- tidy_top(after_p$rotation,  "After")
  
  # Shared symmetric limits across before/after for comparability
  lim <- max(c(before_top$abs_loading, after_top$abs_loading), na.rm = TRUE)
  y_limits <- c(-lim, lim)
  
  # Theme
  big_font_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title   = ggplot2::element_text(size = 12, face = "bold"),
      axis.text    = ggplot2::element_text(size = 9),
      strip.text   = ggplot2::element_text(size = 11, face = "bold"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text  = ggplot2::element_text(size = 10)
    )
  
  # One plotting helper
  mk_plot <- function(df, title) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$variable_wrapped,
                                     y = .data$loading,
                                     fill = .data$sign)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::facet_wrap(~ PC, nrow = 2, scales = "free_y") +
      #ggplot2::scale_y_continuous(limits = y_limits) +
      scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
      ggplot2::scale_fill_manual(values = c("Positive" = "#4C9F50", "Negative" = "#BE0032"),
                                 drop = FALSE) +
      ggplot2::labs(title = title, x = NULL, y = "Loading") +
      big_font_theme +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "none",
        plot.margin = ggplot2::margin(10, 5, 10, 5)
      )
  }
  
  p_before <- mk_plot(before_top, "Before")
  p_after  <- mk_plot(after_top,  "After")
  
  # Standalone legend
  p_leg <- ggplot2::ggplot(before_top, ggplot2::aes(x = variable_wrapped, y = loading, fill = sign)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("Positive" = "#4C9F50", "Negative" = "#BE0032"),
                               name = "Sign", drop = FALSE) +
    ggplot2::guides(y = "none", x = "none", fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    big_font_theme + ggplot2::theme(legend.position = "right")
  leg <- cowplot::get_legend(p_leg)
  
  comb <- cowplot::plot_grid(p_before, p_after, leg,
                             nrow = 1, rel_widths = c(0.7, 0.7, 0.22),
                             labels = NULL, align = "hv", axis = "tblr")
  
  cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("Top ", top_n, " Loadings for PC1 and PC2, Before vs After ", compared_to),
      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5, vjust = 1, size = 14
    ) +
    cowplot::draw_plot(comb, x = 0, y = 0, width = 1, height = 0.93)
}
