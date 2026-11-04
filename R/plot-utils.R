#' Plot helpers, labels, and utilities
#' @keywords internal
#' @noRd
facet_label_map <- function(df) {
  split_df <- split(df, df[["Type"]], drop = TRUE)

  labs <- lapply(split_df, function(x) {
    unique(as.character(x[["Type"]]))[[1L]]
  })

  unlist(labs, use.names = TRUE)
}

#' @keywords internal
#' @noRd
.rsd_axis_range <- function(values) {
  finite_values <- values[is.finite(values)]

  if (!length(finite_values)) {
    return(c(0, 1))
  }

  axis_range <- range(finite_values)
  axis_span <- diff(axis_range)

  if (!is.finite(axis_span) || axis_span == 0) {
    pad <- max(abs(axis_range[[1L]]), 1) * 0.5
    axis_range <- c(axis_range[[1L]] - pad, axis_range[[2L]] + pad)
  }

  axis_range
}

#' @keywords internal
#' @noRd
.rsd_key_label <- function(change, percent) {
  display_change <- dplyr::recode(
    as.character(change),
    "No Change" = "No change",
    .default = as.character(change)
  )

  paste0(display_change, " ", percent, "%")
}

#' @keywords internal
#' @noRd
.rsd_key_rows <- function(d_all, x, y) {
  if (!nrow(d_all)) {
    return(data.frame())
  }

  x_range <- .rsd_axis_range(d_all[[x]])
  y_range <- .rsd_axis_range(d_all[[y]])
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  type_levels <- if (is.factor(d_all[["Type"]])) {
    levels(d_all[["Type"]])
  } else {
    unique(as.character(d_all[["Type"]]))
  }
  split_df <- split(d_all, d_all[["Type"]], drop = TRUE)

  key_x <- x_range[[1L]] + x_span * c(0.005, 0.41, 0.76)
  label_x <- x_range[[1L]] + x_span * c(0.025, 0.435, 0.785)
  label_hjust <- c(0, 0, 0)
  key_y <- y_range[[2L]] + y_span * 0.13
  band_y_min <- y_range[[2L]] + y_span * 0.04
  band_y_max <- y_range[[2L]] + y_span * 0.24
  limit_x <- x_range[[2L]] + x_span * 0.12

  dplyr::bind_rows(lapply(split_df, function(panel_df) {
    panel_type <- unique(as.character(panel_df[["Type"]]))[[1L]]

    pct_tbl(panel_df) |>
      dplyr::mutate(
        Type = factor(panel_type, levels = type_levels),
        label = .rsd_key_label(.data$change, .data$percent),
        key_x = key_x,
        label_x = label_x,
        label_hjust = label_hjust,
        key_y = key_y,
        band_y_min = band_y_min,
        band_y_max = band_y_max,
        limit_x = limit_x,
        limit_y = band_y_max
      ) |>
      dplyr::select(dplyr::all_of(c(
        "Type",
        "change",
        "label",
        "key_x",
        "label_x",
        "label_hjust",
        "key_y",
        "band_y_min",
        "band_y_max",
        "limit_x",
        "limit_y"
      )))
  }))
}

#' @keywords internal
#' @noRd
.rsd_key_band_rows <- function(key_rows) {
  key_rows |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c("Type", "band_y_min", "band_y_max"))))
}

#' @keywords internal
#' @noRd
.rsd_key_limit_rows <- function(key_rows) {
  key_rows |>
    dplyr::transmute(Type = .data$Type, key_x = .data$limit_x, key_y = .data$limit_y) |>
    dplyr::distinct()
}

#' @keywords internal
#' @noRd
mk_plot <- function(d_all, x, y, facet_labs, compared_to) {
  if (!nrow(d_all)) {
    return(ggplot2::ggplot() +
      ggplot2::labs(title = "No points"))
  }

  key_rows <- .rsd_key_rows(d_all, x, y)
  key_band_rows <- .rsd_key_band_rows(key_rows)
  key_limit_rows <- .rsd_key_limit_rows(key_rows)

  ggplot2::ggplot(d_all, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data$change)) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::geom_blank(
      data = key_limit_rows,
      ggplot2::aes(x = .data$key_x, y = .data$key_y),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      data = key_band_rows,
      ggplot2::aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = .data$band_y_min,
        ymax = .data$band_y_max
      ),
      inherit.aes = FALSE,
      fill = "white",
      color = NA
    ) +
    ggplot2::geom_point(
      data = key_rows,
      ggplot2::aes(x = .data$key_x, y = .data$key_y, color = .data$change),
      inherit.aes = FALSE,
      shape = 19,
      size = 2,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = key_rows,
      ggplot2::aes(
        x = .data$label_x,
        y = .data$key_y,
        label = .data$label,
        hjust = .data$label_hjust
      ),
      inherit.aes = FALSE,
      vjust = 0.5,
      size = 3,
      color = "gray20",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = lab_levels,
      labels = c("Increased", "No change", "Decreased"),
      name = "RSD Change"
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0))) +
    ggplot2::facet_wrap(
      ~Type,
      nrow = 1,
      labeller = ggplot2::as_labeller(facet_labs, default = identity),
      strip.position = "top"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey30"),
      strip.text.x = ggplot2::element_text(
        size = 9,
        face = "bold",
        margin = ggplot2::margin(
          t = 6,
          r = 6,
          b = 6,
          l = 6
        )
      ),
      plot.title = ggplot2::element_text(
        size = 14,
        hjust = 0.5,
        face = "bold"
      ),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    ) +
    ggplot2::labs(
      title = paste("Comparison of RSD Before and After", compared_to),
      x = "RSD (%) Before",
      y = "RSD (%) After"
    )
}
pct_tbl <- function(d) {
  total <- nrow(d)
  if (!total) {
    return(stats::setNames(
      data.frame(
        change = factor(lab_levels, levels = lab_levels),
        percent = c(0, 0, 0)
      ),
      c("change", "percent")
    ))
  }
  d |>
    dplyr::count(change, .drop = FALSE) |>
    tidyr::complete(
      change = factor(lab_levels, levels = lab_levels),
      fill = list(n = 0)
    ) |>
    dplyr::mutate(percent = round(100 * n / total, 1)) |>
    dplyr::select(change, percent)
}

#' @keywords internal
#' @noRd
.scatter_metadata_cols <- function() {
  c("sample", "batch", "class", "order")
}

#' @keywords internal
#' @noRd
.scatter_select_metabolite <- function(df, metab) {
  cols <- intersect(c(.scatter_metadata_cols(), metab), names(df))
  df[, cols, drop = FALSE]
}

#' @keywords internal
#' @noRd
.scatter_panel_df <- function(df, panel) {
  df |>
    dplyr::mutate(
      type = ifelse(class == "QC", "QC", "Sample"),
      panel = factor(panel, levels = c("Raw", "Corrected")),
      order = suppressWarnings(as.numeric(order))
    )
}

#' @keywords internal
#' @noRd
.scatter_batch_ranges <- function(df, panel) {
  rng <- df |>
    dplyr::filter(!is.na(order)) |>
    dplyr::group_by(batch) |>
    dplyr::summarise(
      xmin = suppressWarnings(min(as.numeric(order), na.rm = TRUE)),
      xmax = suppressWarnings(max(as.numeric(order), na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(xmin), is.finite(xmax), xmax >= xmin)

  if (nrow(rng) < 2L) {
    return(rng[0, , drop = FALSE])
  }

  rng |>
    dplyr::arrange(xmin) |>
    dplyr::mutate(
      fill = rep(c("lightgray", "white"), length.out = dplyr::n()),
      panel = factor(panel, levels = c("Raw", "Corrected"))
    )
}

#' @keywords internal
#' @noRd
.scatter_prepare_context <- function(data_raw, data_cor, metab) {
  data_raw <- .scatter_select_metabolite(data_raw, metab)
  data_cor <- .scatter_select_metabolite(data_cor, metab)

  raw_panel <- .scatter_panel_df(data_raw, "Raw")
  cor_panel <- .scatter_panel_df(data_cor, "Corrected")

  list(
    data_raw = raw_panel,
    data_cor = cor_panel,
    df_all = dplyr::bind_rows(raw_panel, cor_panel),
    batch_ranges = dplyr::bind_rows(
      .scatter_batch_ranges(raw_panel, "Raw"),
      .scatter_batch_ranges(cor_panel, "Corrected")
    )
  )
}

#' @keywords internal
#' @noRd
.scatter_color_scale <- function() {
  ggplot2::scale_color_manual(
    name = "Type:",
    values = c(Sample = "#F5C710", QC = "#305CDE")
  )
}

#' @keywords internal
#' @noRd
.scatter_base_plot <- function(scatter_context, metab) {
  p <- ggplot2::ggplot(
    scatter_context$df_all,
    ggplot2::aes(x = order, y = .data[[metab]])
  )

  if (nrow(scatter_context$batch_ranges) > 0L) {
    p <- p +
      ggplot2::geom_rect(
        data = scatter_context$batch_ranges,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax, fill = fill),
        ymin = -Inf,
        ymax = Inf,
        inherit.aes = FALSE,
        alpha = 0.3,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_identity(guide = "none")
  }

  p +
    ggplot2::geom_point(
      data = dplyr::filter(scatter_context$df_all, type == "Sample"),
      ggplot2::aes(color = type),
      size = 2,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(scatter_context$df_all, type == "QC"),
      ggplot2::aes(color = type),
      size = 2,
      na.rm = TRUE
    ) +
    .scatter_color_scale() +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(title = metab, x = "Injection Order", y = "Intensity")
}

#' @keywords internal
#' @noRd
.scatter_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.key.width = grid::unit(0.5, "cm"),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      legend.box.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0.5, 0.5, 0.8, 0.5), "cm")
    )
}
