#' Metabolite scatter plot for random forest corrected data
#' @keywords internal
#' @noRd

met_scatter_rf <- function(data_raw, data_cor, i, scatter_context = NULL) {
  scatter_context <- scatter_context %||%
    .scatter_prepare_context(data_raw, data_cor, i)

  # QC SD bands (guard for empty/NA)
  get_stats <- function(df, panel) {
    qc <- df[df$type == "QC", i, drop = TRUE]
    m <- mean(qc, na.rm = TRUE)
    s <- stats::sd(qc, na.rm = TRUE)
    if (!is.finite(m) || !is.finite(s)) {
      return(tibble::tibble(y = numeric(0), sd = factor(), panel = factor()))
    }
    tibble::tibble(
      y     = c(m + s, m - s, m + 2 * s, m - 2 * s),
      sd    = factor(c("\u00B11 SD", "\u00B11 SD", "\u00B12 SD", "\u00B12 SD"), levels = c("\u00B11 SD", "\u00B12 SD")),
      panel = factor(panel, levels = c("Raw", "Corrected"))
    )
  }
  sd_df <- dplyr::bind_rows(
    get_stats(scatter_context$data_raw, "Raw"),
    get_stats(scatter_context$data_cor, "Corrected")
  )

  lty_scale <- ggplot2::scale_linetype_manual(
    name = "SD Range:",
    values = c("\u00B11 SD" = "dashed", "\u00B12 SD" = "solid"),
    guide = ggplot2::guide_legend(override.aes = list(color = c("grey20", "#950606")))
  )

  p <- .scatter_base_plot(scatter_context, i)

  # SD lines only if we have rows
  if (nrow(sd_df)) {
    p <- p +
      ggplot2::geom_hline(
        data = sd_df |> dplyr::filter(sd == "\u00B11 SD"),
        ggplot2::aes(yintercept = y, linetype = sd),
        color = "grey20", linewidth = 0.75
      ) +
      ggplot2::geom_hline(
        data = sd_df |> dplyr::filter(sd == "\u00B12 SD"),
        ggplot2::aes(yintercept = y, linetype = sd),
        color = "#950606", linewidth = 0.75
      )
  }

  p + lty_scale + .scatter_theme()
}
