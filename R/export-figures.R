#' exports figures and returns file path for zip folder
#'
#' @keywords internal
#' @noRd
export_figures <- function(p, d, out_dir = tempdir()) {
  .require_pkg("ggplot2", "write figures")
  fig_dir <- file.path(out_dir, "figures")
  if (dir.exists(fig_dir))
    unlink(fig_dir, recursive = TRUE)
  dir.create(fig_dir)
  
  fmt <- match.arg(p$fig_format, c("png", "pdf"))
  
  mk <- function(subdir) {
    x <- file.path(fig_dir, subdir)
    if (!dir.exists(x))
      dir.create(x, recursive = TRUE)
    x
  }
  
  rsd_dir <- mk("RSD figures")
  pca_dir <- mk("PCA plots")
  met_dir <- mk("metabolite figures")
  
  save_plot <- function(path, plot, w, h) {
    if (fmt == "png") {
      ggplot2::ggsave(
        path,
        plot = plot,
        width = w,
        height = h,
        units = "in",
        dpi = 300,
        bg = "white"
      )
    } else {
      ggplot2::ggsave(
        filename = path,
        plot = plot,
        width = w,
        height = h,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  
  raw_cols <- setdiff(names(d$filtered$df), c("sample", "batch", "class", "order"))
  
  if (isTRUE(p$remove_imputed)) {
    df_cor_mets <- d$filtered_corrected$df_mv
  } else {
    df_cor_mets <- d$filtered_corrected$df_no_mv
  }
  
  cor_cols <- setdiff(names(df_cor_mets), c("sample", "batch", "class", "order"))
  cols <- intersect(raw_cols, cor_cols)
  
  met_paths <- character(0)
  rsd_paths <- character(0)
  pca_paths <- character(0)
  pca_loading_paths <- character(0)
  
  n <- length(cols)
  N <- n + 16
  
  shiny::withProgress(message = "Creating figures...", value = 0, {
    
    rsd_res <- make_all_rsd_plots(p, d)
    for (i in seq_along(rsd_res$rsd_plots)) {
      rsd_path <- file.path(rsd_dir, sprintf("%s.%s", rsd_res$plot_names[i], fmt))
      rsd_paths <- c(rsd_paths, save_plot(rsd_path, rsd_res$rsd_plots[[i]], 7.5, 4.5))
      shiny::incProgress(1 / N, detail = "Saved: rsd figures")
    }
    
    pca_res <- make_all_pca_plots(p, d)
    for (i in seq_along(pca_res$pca_plots)) {
      pca_path <- file.path(pca_dir, sprintf("%s.%s", pca_res$plot_names[i], fmt))
      pca_paths <- c(pca_paths, save_plot(pca_path, pca_res$pca_plots[[i]], 8.333, 4.417))
      shiny::incProgress(1 / N, detail = "Saved: pca figure")
    }
    
    for (i in seq_along(pca_res$pca_loading_plots)) {
      loading_path <- file.path(pca_dir, sprintf("%s.%s", pca_res$loading_plot_names[i], fmt))
      pca_loading_paths <- c(
        pca_loading_paths,
        save_plot(loading_path, pca_res$pca_loading_plots[[i]], 8.75, 4.417)
      )
      shiny::incProgress(1 / N, detail = "Saved: pca loading figure")
    }
    # Export PCA loadings workbook into the PCA plots folder
    pca_loadings_xlsx <- export_pca_loadings_xlsx(
      p = p,
      d = d,
      pca_dir = pca_dir
    )
    shiny::incProgress(1 / N, detail = "Saved: PCA loadings workbook")
    
    for (i in seq_len(n)) {
      metab <- cols[i]
      fig <- make_met_scatter(d, p, metab)
      safe <- gsub("[^A-Za-z0-9_\\-]+", "_", metab)
      path <- file.path(met_dir, sprintf("%s.%s", safe, fmt))
      met_paths <- c(met_paths, save_plot(path, fig, 5, 5))
      shiny::incProgress(1 / n, detail = paste("Saved:", safe))
    }
  })
  
  list(
    fig_dir = normalizePath(fig_dir, winslash = "/", mustWork = TRUE),
    rsd = rsd_paths,
    pca = pca_paths,
    pca_loadings_plots = pca_loading_paths,
    pca_loadings_xlsx = pca_loadings_xlsx,
    metabolite = met_paths
  )
}