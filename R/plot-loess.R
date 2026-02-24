#' metabolite scatter plots for polynomial fit
#' @keywords internal
#' @noRd

met_scatter_loess <- function(data_raw, data_cor, method, i) {
  data_raw <- dplyr::mutate(data_raw, type = ifelse(class=="QC","QC","Sample"),
                            panel = factor("Raw", levels=c("Raw","Corrected")))
  data_cor <- dplyr::mutate(data_cor, type = ifelse(class=="QC","QC","Sample"),
                            panel = factor("Corrected", levels=c("Raw","Corrected")))
  df_all   <- dplyr::bind_rows(data_raw, data_cor)
  
  # ensure order is numeric
  df_all$order    <- suppressWarnings(as.numeric(df_all$order))
  data_raw$order  <- suppressWarnings(as.numeric(data_raw$order))
  data_cor$order  <- suppressWarnings(as.numeric(data_cor$order))
  
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
    
    if (nrow(rng) < 2L) return(rng[0, ])  # skip shading for 0/1 batch
    
    rng |>
      dplyr::arrange(xmin) |>
      dplyr::mutate(
        fill  = rep(c("lightgray","white"), length.out = dplyr::n()),
        panel = factor(panel, levels = c("Raw","Corrected"))
      )
  }
  
  batch_ranges <- dplyr::bind_rows(
    get_batches(data_raw, "Raw"), 
    get_batches(data_cor, "Corrected")
  )
  
  color_scale <- ggplot2::scale_color_manual(
    name = "Type:", values = c(Sample = "#F5C710", QC = "#305CDE")
  )
  
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = order, y = .data[[i]]))
  
  # add background shading only if we have ranges
  if (nrow(batch_ranges)) {
    p <- p + ggplot2::geom_rect(
      data = batch_ranges,
      mapping = ggplot2::aes(xmin = xmin, xmax = xmax, fill = fill),
      ymin = -Inf, ymax = Inf, inherit.aes = FALSE,
      alpha = 0.3, show.legend = FALSE
    ) + ggplot2::scale_fill_identity(guide = "none")
  }
  
  p <- p +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "Sample"),
      ggplot2::aes(color = type), size = 2, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "QC"),
      ggplot2::aes(color = type), size = 2, na.rm = TRUE
    ) +
    color_scale +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(title = i, x = "Injection Order", y = "Intensity") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title   = ggplot2::element_text(size = 12),
      axis.text    = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text  = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      legend.box        = "horizontal",
      legend.box.just   = "center",
      legend.key.width  = grid::unit(0.5, "cm"),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.margin     = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      legend.box.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0.5, 0.5, 0.8, 0.5), "cm")
    )
  
  # LOESS only if enough finite QC points and varying order
  qc_raw <- df_all |>
    dplyr::filter(type == "QC", panel == "Raw") |>
    dplyr::filter(is.finite(order), is.finite(.data[[i]]))
  
  qc_cor <- df_all |>
    dplyr::filter(type == "QC", panel == "Corrected") |>
    dplyr::filter(is.finite(order), is.finite(.data[[i]]))
  
  add_qc_trend <- function(p, df_qc, method, i, span = 0.75) {
    # Need enough points + varying x to draw anything meaningful
    has_line <- nrow(df_qc) >= 3L && dplyr::n_distinct(df_qc$order) >= 2L
    has_ribbon <- nrow(df_qc) >= 10L && dplyr::n_distinct(df_qc$order) >= 3L
    if (!has_line) return(p)
    
    deg <- switch(
      method,
      "local constant regression"    = 0L,
      "local linear regression"      = 1L,  # local linear
      "local polynomial regression" = 2L,  # local polynomial (quadratic)
      2L
    )
    
    print(deg)
    
    p + ggplot2::geom_smooth(
      data = df_qc,
      mapping = ggplot2::aes(order, .data[[i]]),
      method = "loess",
      formula = y ~ x,
      span = span,
      method.args = list(
        degree = deg,
        family = "symmetric",
        control = stats::loess.control(surface = "direct")
      ),
      se = has_ribbon,
      fill = "#305CDE",
      linewidth = 0.75,
      show.legend = FALSE
    )
    
  }
  
  p <- add_qc_trend(p, qc_raw, method = method, i = i, span = 0.75)
  #p <- add_loess(p, qc_cor)
  
  p
}
