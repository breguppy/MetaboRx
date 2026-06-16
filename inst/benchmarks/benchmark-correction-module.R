# Local correction-module benchmark.
# Run with: Rscript inst/benchmarks/benchmark-correction-module.R

set.seed(2)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this benchmark.")
}

pkgload::load_all(helpers = FALSE, quiet = TRUE)

make_correction_data <- function(n_samples = 240L, n_metabolites = 400L) {
  classes <- rep(c("QC", "sample", "sample", "sample"), length.out = n_samples)
  classes[[1L]] <- "QC"
  classes[[n_samples]] <- "QC"

  metab <- matrix(
    stats::rlnorm(n_samples * n_metabolites, meanlog = 10, sdlog = 0.25),
    nrow = n_samples,
    ncol = n_metabolites
  )
  metab[sample(length(metab), size = floor(length(metab) * 0.01))] <- NA_real_
  metab[sample(length(metab), size = floor(length(metab) * 0.005))] <- 0

  df <- data.frame(
    sample = paste0("S", seq_len(n_samples)),
    batch = rep("B1", n_samples),
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

raw_df <- make_correction_data()
metab_cols <- setdiff(names(raw_df), c("sample", "batch", "class", "order"))

imputed <- time_step(
  "impute_missing median/mean",
  impute_missing(raw_df, metab_cols, "median", "mean")
)

lc_corrected <- time_step(
  "correct_data LC",
  correct_data(imputed$df, metab_cols, "LC")
)

time_step(
  "correct_data LOESS",
  correct_data(imputed$df, metab_cols, "LOESS")
)

filtered_corrected <- time_step(
  "filter_by_qc_rsd",
  filter_by_qc_rsd(
    raw_df = raw_df,
    corrected_df = lc_corrected$df,
    rsd_cutoff = 30,
    remove_imputed = FALSE
  )
)

time_step(
  "metabolite_rsd",
  metabolite_rsd(filtered_corrected$df_no_mv)
)

time_step(
  "class_metabolite_rsd",
  class_metabolite_rsd(filtered_corrected$df_no_mv)
)

time_step(
  "transform_data TRN",
  transform_data(
    filtered_corrected = filtered_corrected,
    transform = "TRN",
    withheld_cols = character(0),
    ex_ISTD = TRUE
  )
)

time_step(
  "detect_hotelling_nonqc_dual_z",
  detect_hotelling_nonqc_dual_z(
    df = filtered_corrected$df_no_mv,
    p = list(qcImputeM = "median", samImputeM = "mean"),
    make_pca_plot = FALSE
  )
)

time_step(
  "compute_pairwise_metabolite_correlations",
  compute_pairwise_metabolite_correlations(
    filtered_corrected$df_no_mv,
    cols = setdiff(
      names(filtered_corrected$df_no_mv),
      c("sample", "batch", "class", "order")
    )
  )
)
