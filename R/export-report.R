#' Writes HTML of quality report
#'
#' @keywords internal
#' @noRd
render_report <- function(
  p,
  d,
  out_dir,
  template = system.file("app", "report_templates", "report.Rmd", package = "MetaboRx"),
  rsd_plot_data = NULL,
  pca_pair = NULL,
  hotelling_res = NULL
) {
  .require_pkg("rmarkdown", "render reports")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  env <- new.env(parent = baseenv())

  shiny::withProgress(message = "Creating quality_report.html...", value = 0, {
    met_candidates <- .get_top_two(p, d)
    increased_qc <- .increased_qc_rsd(d)

    met1_plot <- make_met_scatter(d, p, met_candidates[1])
    met2_plot <- make_met_scatter(d, p, met_candidates[2])
    rsd_plot_data <- rsd_plot_data %||% .get_rsd_plot_data(p, d)
    rsd_plot <- make_rsd_plot(
      p,
      d,
      rsd_results = rsd_plot_data$rsd_results,
      compared_to = rsd_plot_data$compared_to
    )

    pca_compare_data <- get_pca_compare_data(
      p = p,
      d = d,
      pca_compare = p$pca_compare %||% "filtered_cor_data"
    )

    meta_df <- NULL
    meta_cols <- c("sample", "batch", "class", "order")

    if (!is.null(d$cleaned) && !is.null(d$cleaned$meta_df)) {
      meta_df <- d$cleaned$meta_df
      meta_cols <- unique(c("sample", setdiff(names(meta_df), "sample")))
    }

    pca_pair <- pca_pair %||%
      compute_pca_pair(
        before = pca_compare_data$before,
        after = pca_compare_data$after,
        p = p,
        before_label = "Before",
        after_label = "After",
        meta_cols = meta_cols,
        meta_df = meta_df,
        sample_col = "sample"
      )

    pca_plot <- plot_pca_from_result(
      p = p,
      pca_pair = pca_pair,
      compared_to = pca_compare_data$compared_to
    )

    pca_loading_plot <- plot_pca_loading_from_result(
      pca_pair = pca_pair,
      compared_to = pca_compare_data$compared_to
    )

    hotelling_res <- hotelling_res %||% d$hotelling_res
    if (is.null(hotelling_res)) {
      hotelling_res <- detect_hotelling_nonqc_dual_z(d$filtered_corrected$df_no_mv, p)
    }
    hotelling_pca_plot <- hotelling_res$pca_plot

    shiny::incProgress(1 / 3, detail = "Saved: plots for report")

    params <- list(
      title = "Metabolomics Data Quality Report",
      notes = p$notes %||% "",
      plots = list(
        "Metabolite Scatter 1" = met1_plot,
        "Metabolite Scatter 2" = met2_plot,
        "RSD Comparison" = rsd_plot,
        "PCA Comparison" = pca_plot,
        "PCA Loading" = pca_loading_plot,
        "Hotelling PCA" = hotelling_pca_plot
      ),
      p = p,
      d = d
    )

    shiny::incProgress(1 / 3, detail = "Saved: report information")

    html_out <- rmarkdown::render(
      input = template,
      output_format = "html_document",
      output_file = "quality_report.html",
      output_dir = out_dir,
      intermediates_dir = out_dir,
      knit_root_dir = out_dir,
      params = params,
      envir = env,
      quiet = TRUE
    )

    shiny::incProgress(1 / 3, detail = "Saved: HTML")

    list(
      html = normalizePath(html_out, winslash = "/")
    )
  })
}
