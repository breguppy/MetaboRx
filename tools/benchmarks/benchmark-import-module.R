# Local import-module benchmark.
# Run with: Rscript tools/benchmarks/benchmark-import-module.R

set.seed(1)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this benchmark.")
}

pkgload::load_all(helpers = FALSE, quiet = TRUE)

make_import_data <- function(n_samples = 600L, n_metabolites = 1500L) {
  classes <- rep(c("QC", "sample", "sample", "sample"), length.out = n_samples)
  classes[[1L]] <- "QC"
  classes[[n_samples]] <- "QC"

  metab <- matrix(
    stats::rlnorm(n_samples * n_metabolites, meanlog = 10, sdlog = 0.25),
    nrow = n_samples,
    ncol = n_metabolites
  )
  metab[sample(length(metab), size = floor(length(metab) * 0.02))] <- NA_real_
  metab[sample(length(metab), size = floor(length(metab) * 0.01))] <- 0

  df <- data.frame(
    sample = paste0("S", seq_len(n_samples)),
    batch = rep(paste0("B", seq_len(6L)), length.out = n_samples),
    class = classes,
    order = seq_len(n_samples),
    check.names = FALSE
  )
  metab_df <- as.data.frame(metab, optional = TRUE)
  names(metab_df) <- paste0("met_", seq_len(n_metabolites))

  cbind(df, metab_df)
}

time_step <- function(label, expr) {
  elapsed <- system.time(result <- force(expr))[["elapsed"]]
  message(sprintf("%s: %.3f sec", label, elapsed))
  result
}

raw_df <- make_import_data()
metab_cols <- setdiff(names(raw_df), c("sample", "batch", "class", "order"))

cleaned <- time_step(
  "clean_data",
  clean_data(
    df = raw_df,
    sample = "sample",
    batch = "batch",
    class = "class",
    order = "order",
    withheld_cols = character(0)
  )
)

time_step(
  "filter_by_missing",
  filter_by_missing(cleaned$df, metab_cols = metab_cols, mv_cutoff = 50)
)

time_step(
  "find_equal_metabolite_cols",
  find_equal_metabolite_cols(cleaned$df, cols = metab_cols, tolerance = 1e-3)
)
