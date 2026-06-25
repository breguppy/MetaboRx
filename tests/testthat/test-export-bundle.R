test_that("export validation names files that are missing or empty", {
  good_file <- tempfile(fileext = ".xlsx")
  empty_file <- tempfile(fileext = ".html")
  missing_file <- tempfile(fileext = ".pdf")
  writeLines("workbook", good_file)
  file.create(empty_file)
  on.exit(unlink(c(good_file, empty_file), force = TRUE), add = TRUE)

  paths <- c(
    "corrected workbook" = good_file,
    "quality report" = empty_file,
    "PCA plot" = missing_file
  )

  expect_error(
    .validate_export_files(paths),
    "quality report, PCA plot",
    fixed = TRUE
  )
})

test_that("export validation accepts nonempty files", {
  files <- c(tempfile(fileext = ".xlsx"), tempfile(fileext = ".pdf"))
  invisible(vapply(files, function(path) {
    writeLines("nonempty", path)
    TRUE
  }, logical(1)))
  on.exit(unlink(files, force = TRUE), add = TRUE)

  expect_invisible(.validate_export_files(files))
})

test_that("figure validation rejects an empty export directory", {
  figure_dir <- tempfile("empty-figures-")
  dir.create(figure_dir)
  on.exit(unlink(figure_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_error(
    .exported_figure_paths(figure_dir),
    "No figure files were generated",
    fixed = TRUE
  )
})

make_complete_export_fixture <- function(fig_format) {
  raw_df <- data.frame(
    sample = paste0("S", seq_len(8)),
    batch = rep(c("B1", "B2"), each = 4),
    class = rep(c("QC", "Case", "Control", "QC"), 2),
    order = seq_len(8),
    met_a = c(10, 15, 18, 12, 11, 16, 20, 13),
    met_b = c(30, 45, 38, 32, 31, 47, 40, 34),
    met_c = c(80, 70, 95, 82, 85, 73, 99, 84),
    check.names = FALSE
  )
  corrected_df <- raw_df
  corrected_df[c("met_a", "met_b", "met_c")] <- lapply(
    corrected_df[c("met_a", "met_b", "met_c")],
    function(values) values / mean(values)
  )

  empty_extremes <- data.frame(
    sample = character(0),
    class = character(0),
    metabolite = character(0),
    z_global = numeric(0),
    abs_z_global = numeric(0),
    z_class = numeric(0),
    abs_z_class = numeric(0),
    T2 = numeric(0),
    check.names = FALSE
  )
  outlier_data <- raw_df[c("sample", "batch", "class", "order")]
  outlier_data$T2 <- seq_len(nrow(outlier_data)) / 10
  outlier_data$is_outlier_sample <- FALSE
  outlier_data$used_in_fit <- outlier_data$class != "QC"
  z_scores <- corrected_df
  z_scores[c("met_a", "met_b", "met_c")] <- 0
  hotelling_res <- list(
    data = outlier_data,
    extreme_values = empty_extremes,
    pc_loadings = data.frame(
      metabolite = c("met_a", "met_b", "met_c"),
      PC1 = c(0.8, 0.4, 0.2),
      PC2 = c(0.1, 0.7, 0.3)
    ),
    z_global = z_scores,
    z_class = z_scores,
    pca_plot = ggplot2::ggplot()
  )

  p <- list(
    sample_col = "sample",
    batch_col = "batch",
    class_col = "class",
    order_col = "order",
    remove_imputed = FALSE,
    rsd_cutoff = 30,
    rsd_filter_threshold = 30,
    post_cor_filter = TRUE,
    remove_qc_average_pct_filter = FALSE,
    transform = "none",
    ex_ISTD = TRUE,
    keep_corrected_qcs = TRUE,
    no_control = TRUE,
    control_class = "",
    withhold_cols = FALSE,
    rsd_compare = "filtered_cor_data",
    rsd_cal = "met",
    rsd_plot_type = "scatter",
    pca_compare = "filtered_cor_data",
    color_col = "class",
    shape_col = "batch",
    fig_format = fig_format,
    qcImputeM = "median",
    samImputeM = "median",
    corr_threshold = 0.8,
    notes = "Automated complete-export test."
  )

  d <- list(
    cleaned = list(
      df = raw_df,
      meta_df = raw_df[c("sample", "batch", "class", "order")],
      replacement_counts = list(
        non_numeric_replaced = 0L,
        zero_replaced = 0L
      ),
      non_numeric_cols = character(0),
      all_missing_zero_qc_cols = character(0),
      duplicate_mets = data.frame(col1 = character(0), col2 = character(0)),
      duplicate_col_names = character(0),
      blank_df = raw_df[0, , drop = FALSE],
      below_blank_threshold_ex_ISTD = character(0),
      withheld_cols = character(0)
    ),
    filtered = list(
      df = raw_df,
      mv_cutoff = 20,
      mv_removed_cols = character(0),
      qc_missing_mets = character(0),
      class_metab_all_missing = character(0),
      blank_threshold = 3,
      remove_blank_threshold_cols = FALSE,
      removed_blank_threshold_cols = character(0),
      blank_threshold_result = NULL
    ),
    imputed = list(
      qc_str = "nothing to impute",
      sam_str = "nothing to impute"
    ),
    corrected = list(
      str = "Random Forest",
      parameters = "the median prediction from three seeded models was used."
    ),
    filtered_corrected = list(
      df_no_mv = corrected_df,
      df_mv = corrected_df,
      rsd_cutoff = 30,
      removed_metabolites_no_mv = character(0),
      removed_metabolites_mv = character(0),
      percent_threshold = 50,
      flagged_mets = character(0),
      removed_mets_pct_diff = character(0)
    ),
    transformed = list(
      df_no_mv = corrected_df,
      df_mv = corrected_df,
      str = "No transformation was applied.",
      withheld_cols = character(0),
      withheld_cols_no_mv = character(0),
      withheld_cols_mv = character(0)
    ),
    hotelling_res = hotelling_res,
    all_corr = NULL
  )

  list(p = p, d = d)
}

if (identical(Sys.getenv("METABORX_RUN_EXPORT_INTEGRATION"), "true")) {
  test_that("complete export bundle contains every real output", {
    expect_true(rmarkdown::pandoc_available())
    fig_format <- Sys.getenv("METABORX_TEST_FIG_FORMAT", unset = "png")
    expect_true(fig_format %in% c("png", "pdf"))
    session <- shiny::MockShinySession$new()
    on.exit(session$close(), add = TRUE)
    fixture <- make_complete_export_fixture(fig_format)
    zip_path <- tempfile(fileext = ".zip")
    on.exit(unlink(zip_path, force = TRUE), add = TRUE)

    result <- suppressWarnings(
      shiny::withReactiveDomain(
        session,
        export_bundle(
          p = fixture$p,
          d = fixture$d,
          file = zip_path,
          export_date = as.Date("2026-01-02")
        )
      )
    )

    expect_true(file.exists(result$file))
    expect_gt(file.size(result$file), 0L)

    contents <- utils::unzip(result$file, list = TRUE)
    archived_files <- contents[!grepl("/$", contents$Name), , drop = FALSE]
    expect_true(all(archived_files$Length > 0L))
    expect_true(any(grepl("^missing_value_counts_.*\\.xlsx$", contents$Name)))
    expect_true(any(grepl("^corrected_data_.*\\.xlsx$", contents$Name)))
    expect_true(any(grepl("^rsd_stats_.*\\.xlsx$", contents$Name)))
    expect_true(any(grepl("^extreme_values_.*\\.xlsx$", contents$Name)))
    expect_true("quality_report.html" %in% contents$Name)
    expect_true(any(grepl("pca_loadings\\.xlsx$", contents$Name)))
    expect_true(any(grepl(paste0("\\.", fig_format, "$"), contents$Name)))
  })
}
