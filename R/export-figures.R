#' Select the most portable available PDF graphics device
#'
#' @keywords internal
#' @noRd
.base_pdf_export_device <- function(filename, ...) {
  grDevices::pdf(file = filename, useDingbats = FALSE, ...)
}

.pdf_export_device <- function(
  cairo_available = getOption(
    "MetaboRx.cairo_available",
    capabilities("cairo")
  )
) {
  if (isTRUE(cairo_available)) {
    grDevices::cairo_pdf
  } else {
    .base_pdf_export_device
  }
}

#' exports figures and returns file path for zip folder
#'
#' @keywords internal
#' @noRd
export_figures <- function(p,
                           d,
                           out_dir = tempdir(),
                           rsd_export = NULL,
                           pca_export = NULL) {
  .require_pkg("ggplot2", "write figures")
  fig_dir <- file.path(out_dir, "figures")
  if (dir.exists(fig_dir)) {
    unlink(fig_dir, recursive = TRUE)
  }
  dir.create(fig_dir)

  fmt <- match.arg(p$fig_format, c("png", "pdf"))

  mk <- function(subdir) {
    x <- file.path(fig_dir, subdir)
    if (!dir.exists(x)) {
      dir.create(x, recursive = TRUE)
    }
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
      pdf_device <- .pdf_export_device()
      save_pdf <- function() {
        ggplot2::ggsave(
          filename = path,
          plot = plot,
          width = w,
          height = h,
          units = "in",
          device = pdf_device
        )
      }

      if (identical(pdf_device, .base_pdf_export_device)) {
        withCallingHandlers(
          save_pdf(),
          warning = function(warning) {
            if (grepl("conversion failure.*mbcsToSbcs", conditionMessage(warning))) {
              invokeRestart("muffleWarning")
            }
          }
        )
      } else {
        save_pdf()
      }
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

  scatter_meta_cols <- intersect(.scatter_metadata_cols(), names(d$filtered$df))
  raw_panel_base <- .scatter_panel_df(d$filtered$df[, scatter_meta_cols, drop = FALSE], "Raw")
  cor_panel_base <- .scatter_panel_df(df_cor_mets[, scatter_meta_cols, drop = FALSE], "Corrected")
  scatter_batch_ranges <- dplyr::bind_rows(
    .scatter_batch_ranges(raw_panel_base, "Raw"),
    .scatter_batch_ranges(cor_panel_base, "Corrected")
  )

  make_export_scatter_context <- function(metab) {
    raw_panel <- dplyr::bind_cols(
      raw_panel_base,
      d$filtered$df[, metab, drop = FALSE]
    )
    cor_panel <- dplyr::bind_cols(
      cor_panel_base,
      df_cor_mets[, metab, drop = FALSE]
    )

    list(
      data_raw = raw_panel,
      data_cor = cor_panel,
      df_all = dplyr::bind_rows(raw_panel, cor_panel),
      batch_ranges = scatter_batch_ranges
    )
  }

  n <- length(cols)
  met_paths <- character(0)
  rsd_paths <- character(0)
  pca_paths <- character(0)
  pca_loading_paths <- character(0)
  pca_loadings_xlsx <- character(0)

  shiny::withProgress(message = "Creating figures...", value = 0, {
    rsd_res <- rsd_export %||% make_all_rsd_plots(p, d)
    pca_res <- pca_export %||%
      suppressWarnings(make_all_pca_plots(p, d, d$cleaned$meta_df))

    total_steps <- length(rsd_res$rsd_plots) +
      length(pca_res$pca_plots) +
      length(pca_res$pca_loading_plots) +
      n +
      1L

    met_paths <- character(n)
    rsd_paths <- character(length(rsd_res$rsd_plots))
    pca_paths <- character(length(pca_res$pca_plots))
    pca_loading_paths <- character(length(pca_res$pca_loading_plots))

    for (i in seq_along(rsd_res$rsd_plots)) {
      rsd_path <- file.path(rsd_dir, sprintf("%s.%s", rsd_res$plot_names[i], fmt))
      rsd_paths[i] <- save_plot(rsd_path, rsd_res$rsd_plots[[i]], 9.375, 5.625)
      shiny::incProgress(1 / total_steps, detail = "Saved: rsd figures")
    }

    for (i in seq_along(pca_res$pca_plots)) {
      pca_path <- file.path(pca_dir, sprintf("%s.%s", pca_res$plot_names[i], fmt))
      pca_paths[i] <- save_plot(pca_path, pca_res$pca_plots[[i]], 12.5, 5.521)
      shiny::incProgress(1 / total_steps, detail = "Saved: pca figure")
    }

    for (i in seq_along(pca_res$pca_loading_plots)) {
      loading_path <- file.path(pca_dir, sprintf("%s.%s", pca_res$loading_plot_names[i], fmt))
      pca_loading_paths[i] <- save_plot(
        loading_path,
        pca_res$pca_loading_plots[[i]],
        10.417,
        5.521
      )
      shiny::incProgress(1 / total_steps, detail = "Saved: pca loading figure")
    }
    # Export PCA loadings workbook into the PCA plots folder
    pca_loadings_xlsx <- export_pca_loadings_xlsx(
      p = p,
      d = d,
      pca_dir = pca_dir,
      pca_results = compute_all_pca_export_results(
        p = p,
        d = d,
        pca_pairs = pca_res$pca_pairs
      )
    )
    shiny::incProgress(1 / total_steps, detail = "Saved: PCA loadings workbook")

    for (i in seq_len(n)) {
      metab <- cols[i]
      suppressWarnings({
        fig <- make_met_scatter(
          d,
          p,
          metab,
          scatter_context = make_export_scatter_context(metab)
        )
        safe <- gsub("[^A-Za-z0-9_\\-]+", "_", metab)
        path <- file.path(met_dir, sprintf("%s.%s", safe, fmt))
        met_paths[i] <- save_plot(path, fig, 6.25, 6.25)
      })
      shiny::incProgress(1 / total_steps, detail = paste("Saved:", safe))
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
