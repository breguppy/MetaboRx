#' Metabolite scatter plot for raw vs corrected data with QC trend overlay
#'
#' Shows raw and corrected metabolite values over injection order. QC points are
#' highlighted and a method-specific trend is overlaid on the raw QC data.
#'
#' For:
#' - local constant regression: plots the actual Nadaraya-Watson fit used for correction
#' - local linear regression: plots a LOESS degree-1 trend
#' - local polynomial regression: plots a LOESS degree-2 trend
#'
#' @param data_raw data.frame
#'   Raw input data.
#' @param data_cor data.frame
#'   Corrected data.
#' @param method character scalar
#'   Correction method description string.
#' @param i character scalar
#'   Metabolite column name to plot.
#' @param span_val numeric scalar
#'   Smoothing span used for trend visualization.
#' @param nw_kernel character scalar
#'   Kernel to use for NW visualization when method is local constant.
#'
#' @return ggplot
#'   Metabolite scatter plot.
#'
#' @keywords internal
#' @noRd
met_scatter_loess <- function(
    data_raw,
    data_cor,
    method,
    i,
    span_val = 0.75,
    nw_kernel = "gaussian"
) {
  data_raw <- dplyr::mutate(
    data_raw,
    type = ifelse(class == "QC", "QC", "Sample"),
    panel = factor("Raw", levels = c("Raw", "Corrected"))
  )
  data_cor <- dplyr::mutate(
    data_cor,
    type = ifelse(class == "QC", "QC", "Sample"),
    panel = factor("Corrected", levels = c("Raw", "Corrected"))
  )
  df_all <- dplyr::bind_rows(data_raw, data_cor)
  
  df_all$order <- suppressWarnings(as.numeric(df_all$order))
  data_raw$order <- suppressWarnings(as.numeric(data_raw$order))
  data_cor$order <- suppressWarnings(as.numeric(data_cor$order))
  
  get_batches <- function(df, panel) {
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
  
  batch_ranges <- dplyr::bind_rows(
    get_batches(data_raw, "Raw"),
    get_batches(data_cor, "Corrected")
  )
  
  color_scale <- ggplot2::scale_color_manual(
    name = "Type:",
    values = c(Sample = "#F5C710", QC = "#305CDE")
  )
  
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = order, y = .data[[i]]))
  
  if (nrow(batch_ranges) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = batch_ranges,
      mapping = ggplot2::aes(xmin = xmin, xmax = xmax, fill = fill),
      ymin = -Inf,
      ymax = Inf,
      inherit.aes = FALSE,
      alpha = 0.3,
      show.legend = FALSE
    ) +
      ggplot2::scale_fill_identity(guide = "none")
  }
  
  p <- p +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "Sample"),
      mapping = ggplot2::aes(color = type),
      size = 2,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "QC"),
      mapping = ggplot2::aes(color = type),
      size = 2,
      na.rm = TRUE
    ) +
    color_scale +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = i,
      x = "Injection Order",
      y = "Intensity"
    ) +
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
  
  qc_raw <- df_all |>
    dplyr::filter(panel == "Raw", type == "QC") |>
    dplyr::filter(is.finite(order), is.finite(.data[[i]]), .data[[i]] > 0)
  
  add_qc_trend <- function(p, df_qc, method, metab, span_val = 0.75, nw_kernel = "gaussian") {
    has_line <- nrow(df_qc) >= 3L && dplyr::n_distinct(df_qc$order) >= 2L
    has_ribbon <- nrow(df_qc) >= 10L && dplyr::n_distinct(df_qc$order) >= 3L
    
    if (!has_line) {
      return(p)
    }
    
    span_val <- as.numeric(span_val)[1]
    if (!is.finite(span_val) || span_val <= 0) {
      span_val <- 0.75
    }
    if (span_val > 1) {
      span_val <- 1
    }
    
    method_key <- tolower(trimws(method))
    
    if (identical(method_key, "local constant regression")) {
      x_grid <- seq(
        min(df_qc$order, na.rm = TRUE),
        max(df_qc$order, na.rm = TRUE),
        length.out = 200L
      )
      
      fit <- .safe_nw_predict_x(
        qc_x = df_qc$order,
        qc_y = df_qc[[metab]],
        newx = x_grid,
        span = span_val,
        kernel = nw_kernel
      )
      
      trend_df <- data.frame(
        order = x_grid,
        fit = fit,
        panel = factor("Raw", levels = c("Raw","Corrected"))
      )
      
      return(
        p + ggplot2::geom_line(
          data = trend_df,
          mapping = ggplot2::aes(x = order, y = fit),
          inherit.aes = FALSE,
          linewidth = 0.75,
          colour = "#305CDE",
          show.legend = FALSE
        )
      )
    }
    
    deg <- switch(
      method_key,
      "local linear regression" = 1L,
      "local polynomial regression" = 2L,
      2L
    )
    
    p + ggplot2::geom_smooth(
      data = df_qc,
      mapping = ggplot2::aes(x = order, y = .data[[metab]]),
      method = "loess",
      formula = y ~ x,
      span = span_val,
      method.args = list(
        degree = deg,
        family = "symmetric"
      ),
      se = has_ribbon,
      fill = "#305CDE",
      colour = "#305CDE",
      linewidth = 0.75,
      show.legend = FALSE
    )
  }
  
  p <- add_qc_trend(
    p = p,
    df_qc = qc_raw,
    method = method,
    metab = i,
    span_val = span_val,
    nw_kernel = nw_kernel
  )
  
  p
}