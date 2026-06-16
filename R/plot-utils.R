#' Plot helpers, labels, and utilities
#' @keywords internal
#' @noRd
facet_label_map <- function(df) {
  split_df <- split(df, df$Type, drop = TRUE)

  labs <- lapply(split_df, function(x) {
    pt <- pct_tbl(x)

    P <- function(k) {
      val <- pt$percent[pt$change == k]
      if (length(val) == 0L) 0 else val
    }

    paste0(
      "<b>", unique(as.character(x$Type)), "</b><br>",
      blk(
        color_values["Increased"],
        "Increased",
        P("Increased"),
        "33.3%",
        "left"
      ),
      blk(
        color_values["No Change"],
        "No change",
        P("No Change"),
        "33.3%",
        "center"
      ),
      blk(
        color_values["Decreased"],
        "Decreased",
        P("Decreased"),
        "33.3%",
        "right"
      )
    )
  })

  unlist(labs, use.names = TRUE)
}

#' @keywords internal
#' @noRd
mk_plot <- function(d_all, x, y, facet_labs, compared_to) {
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop("Install 'ggtext' to use render plots correctly.", call. = FALSE)
  }
  if (!nrow(d_all)) {
    return(ggplot2::ggplot() +
      ggplot2::labs(title = "No points"))
  }
  ggplot2::ggplot(d_all, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = change)) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = lab_levels,
      labels = c("Increased", "No change", "Decreased"),
      name   = "RSD Change"
    ) +
    ggplot2::facet_wrap(
      ~Type,
      nrow = 1,
      labeller = ggplot2::as_labeller(facet_labs, default = identity),
      strip.position = "top"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey30"),
      strip.text.x = ggtext::element_markdown(
        size = 9,
        margin = ggplot2::margin(
          t = 6,
          r = 6,
          b = 6,
          l = 6
        ),
        lineheight = 1.05
      ),
      plot.title = ggplot2::element_text(
        size = 14,
        hjust = 0.5,
        face = "bold"
      ),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
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

#' @keywords internal
#' @noRd
blk <- function(col,
                lab,
                pct,
                width = "33.3%",
                align = "left") {
  sprintf(
    "<span style='display:inline-block; width:%s; text-align:%s; white-space:nowrap'>
            %s&nbsp;%s&nbsp;%s%%
          </span>",
    width,
    align,
    circle(col),
    lab,
    pct
  )
}

#' @keywords internal
#' @noRd
circle <- function(col, ptsize = 12) {
  sprintf(
    "<span style='color:%s; font-size:%dpt'>&#9679;</span>",
    col,
    ptsize
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
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
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
