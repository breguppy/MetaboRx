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
  nw_kernel = "gaussian",
  scatter_context = NULL
) {
  scatter_context <- scatter_context %||%
    .scatter_prepare_context(data_raw, data_cor, i)

  p <- .scatter_base_plot(scatter_context, i) +
    .scatter_theme()

  qc_raw <- scatter_context$df_all |>
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
        panel = factor("Raw", levels = c("Raw", "Corrected"))
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

    deg <- switch(method_key,
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
