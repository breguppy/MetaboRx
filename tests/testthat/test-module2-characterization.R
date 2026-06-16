test_that("percent distance filtering returns flagged and retained metabolites", {
  df <- data.frame(
    sample = paste0("s", 1:6),
    batch = 1L,
    class = c("QC", "QC", "QC", "sample", "sample", "sample"),
    order = 1:6,
    near_qc = c(10, 10, 10, 11, 10, 9),
    far_qc = c(10, 10, 10, 30, 30, 30),
    zero_qc = c(0, 0, 0, 1, 1, 1),
    stringsAsFactors = FALSE
  )

  stats <- get_metabs_pct_diff_vs_qc_average(
    df,
    percent_threshold = 100,
    return_stats = TRUE
  )
  removed <- remove_metabs_pct_diff_vs_qc_average(
    df,
    percent_threshold = 100,
    return_result = TRUE
  )

  expect_named(stats, c(
    "metabolite",
    "sample_mean",
    "qc_mean",
    "percent_distance_from_qc_average",
    "flagged"
  ))
  expect_setequal(stats$metabolite[stats$flagged], c("far_qc", "zero_qc"))
  expect_setequal(removed$removed_metabolites, c("far_qc", "zero_qc"))
  expect_setequal(names(removed$df), c("sample", "batch", "class", "order", "near_qc"))
})

test_that("module 2 helper chain preserves expected list contracts", {
  df <- data.frame(
    sample = paste0("s", 1:7),
    batch = "A",
    class = c("QC", "sample", "QC", "sample", "QC", "sample", "QC"),
    order = 1:7,
    M1 = c(100, 105, 110, 115, 120, 125, 130),
    M2 = c(50, NA, 55, 60, 65, 70, 75),
    ISTD_M1 = c(20, 21, 20, 22, 21, 22, 21),
    stringsAsFactors = FALSE
  )
  metab_cols <- c("M1", "M2", "ISTD_M1")

  imputed <- impute_missing(df, metab_cols, "median", "mean")
  corrected <- correct_data(imputed$df, metab_cols, "LC")
  filtered <- filter_by_qc_rsd(
    raw_df = df,
    corrected_df = corrected$df,
    rsd_cutoff = Inf,
    remove_imputed = FALSE
  )
  transformed <- transform_data(
    filtered_corrected = filtered,
    transform = "none",
    withheld_cols = character(0),
    ex_ISTD = TRUE
  )

  expect_named(imputed, c("df", "qc_str", "sam_str", "n_missv"))
  expect_named(corrected, c("df", "str", "parameters"))
  expect_named(filtered, c(
    "df_no_mv",
    "df_mv",
    "rsd_cutoff",
    "removed_metabolites_no_mv",
    "removed_metabolites_mv"
  ))
  expect_named(transformed, c(
    "df_mv",
    "df_no_mv",
    "equal_weight_df_mv",
    "equal_weight_df_no_mv",
    "str",
    "withheld_cols_mv",
    "withheld_cols_no_mv"
  ))
  expect_setequal(transformed$withheld_cols_no_mv, "ISTD_M1")
})
