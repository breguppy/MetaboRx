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
