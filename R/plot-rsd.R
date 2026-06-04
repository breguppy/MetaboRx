# ------------------------------------------------------------------------------
# Scatter comparison plots
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
plot_rsd_comparison <- function(df_before, df_after, compared_to) {
  compare_df <- .build_rsd_results(df_before, df_after)$metabolite$compare
  
  d_all <- compare_df |>
    dplyr::filter(is.finite(before), is.finite(after)) |>
    dplyr::select(Type, before, after, change)
  
  if (!nrow(d_all)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No overlapping metabolites"))
  }
  
  facet_labs <- facet_label_map(d_all)
  mk_plot(d_all, "before", "after", facet_labs, compared_to)
}

#' @keywords internal
#' @noRd
plot_rsd_comparison_class_met <- function(df_before, df_after, compared_to) {
  compare_df <- .build_rsd_results(df_before, df_after)$class_metabolite$compare
  
  d_all <- compare_df |>
    dplyr::filter(is.finite(before), is.finite(after)) |>
    dplyr::select(Type, before, after, change)
  
  if (!nrow(d_all)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No overlapping metabolites"))
  }
  
  facet_labs <- facet_label_map(d_all)
  mk_plot(d_all, "before", "after", facet_labs, compared_to)
}

# ------------------------------------------------------------------------------
# Density plots
# ------------------------------------------------------------------------------

#' Plot metabolite RSD distributions before and after correction
#'
#' @param df_before Data frame for computing raw RSD.
#' @param df_after Data frame for computing metabolite-level RSD after correction
#'   or correction and transformation.
#' @param compared_to Character. Type of comparison, usually "Correction" or
#'   "Correction and Transformation".
#' @param before_label Character. Label for the before group in the legend.
#' @param after_label Character. Label for the after group in the legend.
#' @param facet_scales Character. Facet scale behavior passed to
#'   `ggplot2::facet_wrap()`. Use `"free_x"` to allow QC and Samples panels to
#'   use separate x-axis ranges. Use `"fixed"` to force the same x-axis range.
#'
#' @return A ggplot object with two panels: Samples and QC.
#' @keywords internal
#' @noRd
plot_met_rsd_distributions <- function(df_before,
                                       df_after,
                                       compared_to,
                                       before_label = "Before",
                                       after_label = "After",
                                       facet_scales = "free_x") {
  valid_scales <- c("fixed", "free", "free_x", "free_y")
  
  if (!facet_scales %in% valid_scales) {
    stop(
      "'facet_scales' must be one of: ",
      paste(valid_scales, collapse = ", "),
      call. = FALSE
    )
  }
  
  rsd_results <- .build_rsd_results(df_before, df_after)
  rsd_before <- rsd_results$metabolite$before
  rsd_after <- rsd_results$metabolite$after
  
  rsd_before2 <- dplyr::mutate(rsd_before, dataset = before_label)
  rsd_after2 <- dplyr::mutate(rsd_after, dataset = after_label)
  
  rsd_all <- dplyr::bind_rows(rsd_before2, rsd_after2)
  
  rsd_long <- tidyr::pivot_longer(
    rsd_all,
    cols = c("RSD_QC", "RSD_NonQC"),
    names_to = "type",
    values_to = "RSD"
  )
  
  rsd_long$type <- factor(
    rsd_long$type,
    levels = c("RSD_NonQC", "RSD_QC"),
    labels = c("Samples", "QC")
  )
  
  rsd_long$dataset <- factor(
    rsd_long$dataset,
    levels = c(before_label, after_label)
  )
  
  rsd_long <- rsd_long[!is.na(rsd_long$RSD), , drop = FALSE]
  
  col_vals <- stats::setNames(
    c("#1F77B4", "#FF7F0E"),
    c(before_label, after_label)
  )
  
  ggplot2::ggplot(
    rsd_long,
    ggplot2::aes(x = RSD, fill = dataset, color = dataset)
  ) +
    ggplot2::facet_wrap(~ type, nrow = 1, scales = facet_scales) +
    ggplot2::geom_density(
      data = subset(rsd_long, dataset == before_label),
      alpha = 0.3,
      adjust = 1
    ) +
    ggplot2::geom_density(
      data = subset(rsd_long, dataset == after_label),
      alpha = 0.3,
      adjust = 1
    ) +
    ggplot2::scale_fill_manual(values = col_vals) +
    ggplot2::scale_color_manual(values = col_vals) +
    ggplot2::labs(
      title = paste("Comparison of RSD Before and After", compared_to),
      x = "RSD (%)",
      y = "Density",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey30"),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    )
}


#' Plot class metabolite RSD distributions before and after correction
#'
#' @param df_before Data frame for computing raw RSD.
#' @param df_after Data frame for computing class-metabolite RSD after correction
#'   or correction and transformation.
#' @param compared_to Character. Type of comparison, usually "Correction" or
#'   "Correction and Transformation".
#' @param before_label Character. Label for the before group in the legend.
#' @param after_label Character. Label for the after group in the legend.
#' @param facet_scales Character. Facet scale behavior passed to
#'   `ggplot2::facet_wrap()`. Use `"free_x"` to allow QC and Samples panels to
#'   use separate x-axis ranges. Use `"fixed"` to force the same x-axis range.
#'
#' @return A ggplot object with two panels: Samples and QC.
#' @keywords internal
#' @noRd
plot_class_rsd_distributions <- function(df_before,
                                         df_after,
                                         compared_to,
                                         before_label = "Before",
                                         after_label = "After",
                                         facet_scales = "free_x") {
  valid_scales <- c("fixed", "free", "free_x", "free_y")
  
  if (!facet_scales %in% valid_scales) {
    stop(
      "'facet_scales' must be one of: ",
      paste(valid_scales, collapse = ", "),
      call. = FALSE
    )
  }
  
  rsd_results <- .build_rsd_results(df_before, df_after)
  rsd_before <- rsd_results$class_metabolite$before
  rsd_after <- rsd_results$class_metabolite$after
  
  rsd_before2 <- dplyr::mutate(rsd_before, dataset = before_label)
  rsd_after2 <- dplyr::mutate(rsd_after, dataset = after_label)
  
  rsd_all <- dplyr::bind_rows(rsd_before2, rsd_after2)
  
  rsd_all$dataset <- factor(
    rsd_all$dataset,
    levels = c(before_label, after_label)
  )
  
  rsd_all$type <- ifelse(rsd_all$class == "QC", "QC", "Samples")
  rsd_all$type <- factor(rsd_all$type, levels = c("Samples", "QC"))
  
  rsd_all <- rsd_all[!is.na(rsd_all$RSD), , drop = FALSE]
  
  col_vals <- stats::setNames(
    c("#1F77B4", "#FF7F0E"),
    c(before_label, after_label)
  )
  
  ggplot2::ggplot(
    rsd_all,
    ggplot2::aes(x = RSD, fill = dataset, color = dataset)
  ) +
    ggplot2::facet_wrap(~ type, nrow = 1, scales = facet_scales) +
    ggplot2::geom_density(
      data = subset(rsd_all, dataset == before_label),
      alpha = 0.3,
      adjust = 1
    ) +
    ggplot2::geom_density(
      data = subset(rsd_all, dataset == after_label),
      alpha = 0.3,
      adjust = 1
    ) +
    ggplot2::scale_fill_manual(values = col_vals) +
    ggplot2::scale_color_manual(values = col_vals) +
    ggplot2::labs(
      title = paste("Comparison of RSD Before and After", compared_to),
      x = "RSD (%)",
      y = "Density",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey30"),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    )
}