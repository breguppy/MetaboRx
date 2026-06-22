#' Validate generated export files before archiving
#'
#' @keywords internal
#' @noRd
.validate_export_files <- function(paths) {
  labels <- names(paths)
  paths <- as.character(paths)
  if (is.null(labels)) {
    labels <- basename(paths)
  } else {
    labels[!nzchar(labels)] <- basename(paths[!nzchar(labels)])
  }

  exists <- file.exists(paths)
  sizes <- rep(NA_real_, length(paths))
  sizes[exists] <- file.size(paths[exists])
  invalid <- !exists | is.na(sizes) | sizes <= 0L

  if (any(invalid)) {
    stop(
      "Export could not be completed. Missing or empty files: ",
      paste(labels[invalid], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(paths)
}

.exported_figure_paths <- function(figure_dir) {
  if (length(figure_dir) != 1L || !dir.exists(figure_dir)) {
    stop("The figure export directory was not created.", call. = FALSE)
  }

  paths <- list.files(figure_dir, recursive = TRUE, full.names = TRUE)
  if (length(paths) == 0L) {
    stop("No figure files were generated.", call. = FALSE)
  }

  .validate_export_files(
    stats::setNames(paths, paste0("figure ", seq_along(paths)))
  )
  paths
}

#' Build and validate the complete MetaboRx download archive
#'
#' @keywords internal
#' @noRd
export_bundle <- function(p, d, file, export_date = Sys.Date()) {
  .require_pkg("zip", "create a zip archive")

  export_cache <- .export_bundle_cache()
  base_dir <- tempfile("bundle_")
  dir.create(base_dir)
  on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)

  mv_xlsx_path <- file.path(
    base_dir,
    sprintf("missing_value_counts_%s.xlsx", export_date)
  )
  mv_wb <- export_mv_xlsx(p, d)
  openxlsx::saveWorkbook(mv_wb, mv_xlsx_path, overwrite = TRUE)

  xlsx_path <- file.path(
    base_dir,
    sprintf("corrected_data_%s.xlsx", export_date)
  )
  wb <- export_xlsx(p, d)
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  temp_params <- p
  temp_params$rsd_compare <- "filtered_cor_data"
  stats_xlsx_path1 <- file.path(
    base_dir,
    sprintf("rsd_stats_%s_%s.xlsx", temp_params$rsd_compare, export_date)
  )
  stats_wb1 <- export_stats_xlsx(
    temp_params,
    d,
    rsd_results = .export_cached_rsd_results(
      export_cache,
      "filtered_cor_data",
      p,
      d
    )
  )
  openxlsx::saveWorkbook(stats_wb1, stats_xlsx_path1, overwrite = TRUE)

  stats_xlsx_path2 <- character(0)
  if (!identical(p$transform, "none")) {
    temp_params$rsd_compare <- "transformed_cor_data"
    stats_xlsx_path2 <- file.path(
      base_dir,
      sprintf("rsd_stats_%s_%s.xlsx", temp_params$rsd_compare, export_date)
    )
    stats_wb2 <- export_stats_xlsx(
      temp_params,
      d,
      rsd_results = .export_cached_rsd_results(
        export_cache,
        "transformed_cor_data",
        p,
        d
      )
    )
    openxlsx::saveWorkbook(stats_wb2, stats_xlsx_path2, overwrite = TRUE)
  }

  corr_xlsx_path <- character(0)
  if (!is.null(d$all_corr)) {
    corr_xlsx_path <- file.path(
      base_dir,
      sprintf("metabolite_correlations_%s.xlsx", export_date)
    )
    corr_wb <- export_corr_xlsx(d$all_corr)
    openxlsx::saveWorkbook(corr_wb, corr_xlsx_path, overwrite = TRUE)
  }

  outlier_xlsx_path <- file.path(
    base_dir,
    sprintf("extreme_values_%s.xlsx", export_date)
  )
  outlier_wb <- export_outliers_xlsx(
    p,
    d,
    hotelling_res = .export_cached_hotelling(export_cache, p, d)
  )
  openxlsx::saveWorkbook(outlier_wb, outlier_xlsx_path, overwrite = TRUE)

  figs <- export_figures(
    p,
    d,
    out_dir = base_dir,
    rsd_export = .export_cached_rsd_plots(export_cache, p, d),
    pca_export = .export_cached_pca_plots(export_cache, p, d)
  )

  render_report(
    p,
    d,
    out_dir = base_dir,
    rsd_plot_data = .export_cached_rsd_plot_data(
      export_cache,
      p$rsd_compare,
      p,
      d
    ),
    pca_pair = .export_cached_pca_pair(export_cache, p, d),
    hotelling_res = .export_cached_hotelling(export_cache, p, d)
  )

  figure_paths <- .exported_figure_paths(figs$fig_dir)
  report_path <- file.path(base_dir, "quality_report.html")
  generated_files <- c(
    "missing-value workbook" = mv_xlsx_path,
    "corrected-data workbook" = xlsx_path,
    "corrected RSD workbook" = stats_xlsx_path1,
    if (length(stats_xlsx_path2)) {
      c("transformed RSD workbook" = stats_xlsx_path2)
    },
    if (length(corr_xlsx_path)) {
      c("correlation workbook" = corr_xlsx_path)
    },
    "extreme-value workbook" = outlier_xlsx_path,
    stats::setNames(figure_paths, paste0("figure ", seq_along(figure_paths))),
    "quality report" = report_path
  )
  .validate_export_files(generated_files)

  relative_entries <- c(
    "figures",
    basename(mv_xlsx_path),
    basename(xlsx_path),
    basename(stats_xlsx_path1),
    basename(outlier_xlsx_path),
    basename(corr_xlsx_path),
    "quality_report.html",
    basename(stats_xlsx_path2)
  )

  temp_zip <- tempfile(fileext = ".zip")
  on.exit(unlink(temp_zip, force = TRUE), add = TRUE)
  zip::zipr(zipfile = temp_zip, files = relative_entries, root = base_dir)
  .validate_export_files(c("ZIP archive" = temp_zip))

  copied <- file.copy(temp_zip, file, overwrite = TRUE)
  if (!isTRUE(copied)) {
    stop("Export archive could not be copied to the download location.", call. = FALSE)
  }
  .validate_export_files(c("download archive" = file))

  invisible(list(
    file = normalizePath(file, winslash = "/", mustWork = TRUE),
    entries = relative_entries,
    generated_files = generated_files
  ))
}
