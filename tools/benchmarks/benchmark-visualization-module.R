# Visualization-module benchmark.
# Run with: Rscript tools/benchmarks/benchmark-visualization-module.R
#
# Suggested workflow:
# 1. Run this script before an optimization change and save the console output.
# 2. Run it again after the change with the same n_samples/n_metabolites values.
# 3. Record before/after elapsed times in the PR notes.

set.seed(3)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this benchmark.")
}

pkgload::load_all(helpers = FALSE, quiet = TRUE)

make_visualization_data <- function(n_samples = 180L, n_metabolites = 250L) {
  classes <- rep(c("QC", "sample", "sample", "sample", "sample"), length.out = n_samples)
  classes[[1L]] <- "QC"
  classes[[n_samples]] <- "QC"

  metab <- matrix(
    stats::rlnorm(n_samples * n_metabolites, meanlog = 10, sdlog = 0.25),
    nrow = n_samples,
    ncol = n_metabolites
  )

  df <- data.frame(
    sample = paste0("S", seq_len(n_samples)),
    batch = rep(paste0("B", seq_len(3)), length.out = n_samples),
    class = classes,
    order = seq_len(n_samples),
    check.names = FALSE
  )

  metab_df <- as.data.frame(metab, optional = TRUE)
  names(metab_df) <- paste0("met_", seq_len(n_metabolites))
  raw_df <- cbind(df, metab_df)

  corrected_df <- raw_df
  corrected_df[, names(metab_df)] <- corrected_df[, names(metab_df)] |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / stats::median(.x, na.rm = TRUE)))

  transformed_df <- corrected_df
  transformed_df[, names(metab_df)] <- transformed_df[, names(metab_df)] |>
    dplyr::mutate(dplyr::across(dplyr::everything(), log1p))

  list(
    filtered = list(df = raw_df),
    corrected = list(str = "Random Forest"),
    filtered_corrected = list(df_no_mv = corrected_df, df_mv = corrected_df),
    transformed = list(df_no_mv = transformed_df, df_mv = transformed_df),
    cleaned = list(meta_df = raw_df[, c("sample", "batch", "class", "order")])
  )
}

time_step <- function(label, expr) {
  elapsed <- system.time(result <- force(expr))[["elapsed"]]
  message(sprintf("%s: %.3f sec", label, elapsed))
  result
}

p <- list(
  remove_imputed = FALSE,
  transform = "log",
  fig_format = "png",
  qcImputeM = "median",
  samImputeM = "median",
  rsd_compare = "filtered_cor_data",
  rsd_cal = "met",
  rsd_plot_type = "dist",
  pca_compare = "filtered_cor_data",
  color_col = "class",
  shape_col = "batch"
)

d <- make_visualization_data()
metab_cols <- setdiff(names(d$filtered$df), c("sample", "batch", "class", "order"))

time_step(
  "make_met_scatter single metabolite",
  make_met_scatter(d, p, metab_cols[[1L]])
)

time_step(
  "all-metabolite scatter context and plot loop",
  lapply(metab_cols, function(metab) {
    make_met_scatter(
      d,
      p,
      metab,
      scatter_context = .scatter_prepare_context(
        d$filtered$df,
        d$filtered_corrected$df_no_mv,
        metab
      )
    )
  })
)

time_step(
  "make_all_rsd_plots",
  make_all_rsd_plots(p, d)
)

time_step(
  "make_all_pca_plots",
  make_all_pca_plots(p, d, d$cleaned$meta_df)
)

if (
  requireNamespace("cowplot", quietly = TRUE) &&
    requireNamespace("ggtext", quietly = TRUE) &&
    requireNamespace("openxlsx", quietly = TRUE)
) {
  out_dir <- tempfile("metaborx-viz-benchmark-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  time_step(
    "export_figures png",
    export_figures(p, d, out_dir = out_dir)
  )
} else {
  message("Skipping export_figures png: cowplot, ggtext, and openxlsx are required.")
}
