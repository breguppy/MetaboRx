make_export_figures_data <- function() {
  df <- data.frame(
    sample = paste0("S", seq_len(6)),
    batch = rep(c("B1", "B2"), each = 3),
    class = c("QC", "Sample", "Sample", "QC", "Sample", "Sample"),
    order = seq_len(6),
    met_a = c(10, 11, 12, 10.5, 11.5, 12.5),
    met_b = c(20, 21, 19, 20.5, 21.5, 19.5),
    check.names = FALSE
  )

  list(
    filtered = list(df = df),
    corrected = list(str = "Random Forest"),
    filtered_corrected = list(df_no_mv = df, df_mv = df),
    transformed = list(df_no_mv = df, df_mv = df),
    cleaned = list(meta_df = df[, c("sample", "batch", "class", "order")])
  )
}

test_that("export_figures preserves returned file groups and names", {
  testthat::skip_if_not_installed("cowplot")
  testthat::skip_if_not_installed("ggtext")
  testthat::skip_if_not_installed("openxlsx")

  out_dir <- tempfile("figures-test-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  p <- list(
    remove_imputed = FALSE,
    transform = "none",
    fig_format = "png",
    qcImputeM = "median",
    samImputeM = "median"
  )

  session <- shiny::MockShinySession$new()
  res <- shiny::withReactiveDomain(
    session,
    export_figures(p = p, d = make_export_figures_data(), out_dir = out_dir)
  )

  testthat::expect_setequal(
    names(res),
    c(
      "fig_dir",
      "rsd",
      "pca",
      "pca_loadings_plots",
      "pca_loadings_xlsx",
      "metabolite"
    )
  )
  testthat::expect_true(dir.exists(res$fig_dir))
  testthat::expect_true(all(file.exists(unlist(res[c("rsd", "pca", "pca_loadings_plots", "metabolite")]))))
  testthat::expect_true(file.exists(res$pca_loadings_xlsx))
  testthat::expect_true(all(grepl("\\.png$", c(res$rsd, res$pca, res$pca_loadings_plots, res$metabolite))))
  testthat::expect_true(grepl("pca_loadings\\.xlsx$", res$pca_loadings_xlsx))
})
