make_export_xlsx_fixture <- function() {
  raw_df <- data.frame(
    sample = c("QC1", "S1", "S2", "QC2"),
    batch = "B1",
    class = c("QC", "Sample", "Sample", "QC"),
    order = 1:4,
    met_a = c(1234.56789, 2.34567, 3.45678, 4.56789),
    met_b = c(1, 2, 3, 100),
    check.names = FALSE
  )

  transformed_df <- raw_df
  transformed_df$met_a <- c(1234.56789, 2.34567, 3.45678, 4.56789)
  transformed_df$met_b <- c(1.23456, 2.34567, 3.45678, 100.98765)

  p <- list(
    sample_col = "sample",
    batch_col = "batch",
    class_col = "class",
    order_col = "order",
    remove_imputed = FALSE,
    rsd_cutoff = Inf,
    rsd_filter_threshold = 30,
    post_cor_filter = TRUE,
    remove_qc_average_pct_filter = FALSE,
    transform = "none",
    ex_ISTD = TRUE,
    keep_corrected_qcs = TRUE,
    no_control = TRUE,
    control_class = "",
    withhold_cols = TRUE
  )

  d <- list(
    cleaned = list(
      meta_df = raw_df[, c("sample", "batch", "class", "order"), drop = FALSE],
      withheld_cols = "operator_note",
      all_missing_zero_qc_cols = "all_zero_qc",
      duplicate_col_names = "duplicate_raw",
      duplicate_mets = data.frame(
        col1 = "dup_a",
        col2 = "dup_b",
        stringsAsFactors = FALSE
      )
    ),
    filtered = list(
      df = raw_df,
      mv_cutoff = 20,
      mv_removed_cols = "missing_removed",
      qc_missing_mets = "qc_missing_after_filter",
      blank_threshold = 3,
      remove_blank_threshold_cols = FALSE,
      removed_blank_threshold_cols = character(0),
      blank_threshold_result = list(
        below_blank_threshold_ex_ISTD = "blank_flagged"
      )
    ),
    imputed = list(
      qc_str = "No QC imputation",
      sam_str = "No sample imputation"
    ),
    corrected = list(
      str = "No correction",
      parameters = "No correction parameters."
    ),
    filtered_corrected = list(
      df_no_mv = transformed_df,
      df_mv = transformed_df,
      rsd_cutoff = Inf,
      removed_metabolites_no_mv = character(0),
      removed_metabolites_mv = character(0),
      percent_threshold = 50,
      flagged_mets = "average_flagged",
      removed_mets_pct_diff = character(0)
    ),
    transformed = list(
      df_no_mv = transformed_df,
      df_mv = transformed_df,
      withheld_cols_no_mv = character(0),
      withheld_cols_mv = character(0)
    )
  )

  list(p = p, d = d)
}

export_xlsx_fixture <- function(fixture, file) {
  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(
    session,
    export_xlsx(fixture$p, fixture$d, file = file)
  )
}

test_that("export_xlsx rounds metabolite values only in non-raw sheets", {
  fixture <- make_export_xlsx_fixture()
  original <- fixture$d$transformed$df_no_mv
  file <- tempfile(fileext = ".xlsx")

  export_xlsx_fixture(fixture, file)

  raw <- openxlsx::read.xlsx(file, sheet = "0. Raw Data", startRow = 3)
  corrected <- openxlsx::read.xlsx(file, sheet = "2. Drift Corrected", startRow = 3)
  normalized <- openxlsx::read.xlsx(file, sheet = "3. Samples Normalized", startRow = 3)
  metaboanalyst <- openxlsx::read.xlsx(file, sheet = "Appendix1. MetaboAnalyst Ready")

  expect_equal(raw$met_a[[1]], 1234.56789)
  expect_equal(corrected$met_a[[1]], 1234.568)
  expect_equal(normalized$met_a[[1]], 1234.568)
  expect_equal(metaboanalyst$met_a[[1]], 1234.568)
  expect_equal(fixture$d$transformed$df_no_mv, original)
})

test_that("export_xlsx writes expanded correction settings audit tables", {
  fixture <- make_export_xlsx_fixture()
  file <- tempfile(fileext = ".xlsx")

  export_xlsx_fixture(fixture, file)

  settings <- openxlsx::read.xlsx(
    file,
    sheet = "1. Correction Settings",
    colNames = FALSE
  )
  settings_values <- as.character(unlist(settings, use.names = FALSE))

  expect_true("Blank Threshold Multiplier" %in% settings_values)
  expect_true("Sample/QC Average Difference Threshold" %in% settings_values)
  expect_true("All Missing/Zero in QC Metabolites Removed" %in% settings_values)
  expect_true("Equal/Duplicate Metabolite Pairs Detected" %in% settings_values)
  expect_true("Blank-Threshold Flagged Metabolites" %in% settings_values)
  expect_true("Sample/QC Average Flagged Metabolites" %in% settings_values)
  expect_true("QC-RSD Flagged Metabolites" %in% settings_values)
  expect_true("QC Missing After Missing-Value Filtering" %in% settings_values)
  expect_true("blank_flagged" %in% settings_values)
  expect_true("average_flagged" %in% settings_values)
  expect_true("met_b" %in% settings_values)
})
