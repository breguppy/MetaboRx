make_cached_export_df <- function() {
  data.frame(
    sample = paste0("S", seq_len(8)),
    batch = rep(c("B1", "B2"), each = 4),
    class = c("QC", "QC", "Case", "Control", "QC", "QC", "Case", "Control"),
    order = seq_len(8),
    met_a = c(10, 12, 30, 35, 11, 13, 31, 37),
    met_b = c(50, 48, 70, 72, 49, 47, 69, 73),
    check.names = FALSE
  )
}

make_cached_hotelling_result <- function() {
  data <- data.frame(
    sample = c("S1", "S2"),
    batch = "B1",
    class = c("Case", "QC"),
    order = 1:2,
    T2 = c(12, 1),
    is_outlier_sample = c(TRUE, FALSE),
    used_in_fit = c(TRUE, FALSE)
  )
  
  list(
    data = data,
    extreme_values = data.frame(
      sample = data$sample,
      batch = data$batch,
      class = data$class,
      order = data$order,
      metabolite = "met_a",
      method = c("Grubbs-type", "ROUT"),
      test_statistic = c(6.1234, -9.5678),
      adjusted_p = c(0.002, 0.0005),
      p_cutoff = c(0.005, 0.001),
      T2 = data$T2
    ),
    params = list(
      min_rout_observations = 8L,
      grubbs_alpha = 0.005,
      rout_q = 0.001,
      log_offset = 1,
      class_col = "class",
      qc_label = "QC"
    )
  )
}


test_that("export_stats_xlsx reuses supplied RSD results", {
  df_before <- make_cached_export_df()
  df_after <- dplyr::mutate(
    df_before,
    met_a = .data$met_a / mean(.data$met_a),
    met_b = .data$met_b / mean(.data$met_b)
  )
  rsd_results <- .build_rsd_results(df_before, df_after)
  file <- tempfile(fileext = ".xlsx")

  testthat::local_mocked_bindings(
    .build_rsd_results = function(...) {
      stop("RSD results should be reused")
    }
  )

  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(
    session,
    export_stats_xlsx(
      list(rsd_compare = "filtered_cor_data", remove_imputed = FALSE),
      list(
        filtered = list(df = df_before),
        filtered_corrected = list(df_no_mv = df_after, df_mv = df_after)
      ),
      file = file,
      rsd_results = rsd_results
    )
  )

  comparison <- openxlsx::read.xlsx(file, sheet = "Metabolite RSD Comparison", startRow = 3)
  expect_setequal(names(comparison), c("Metabolite", "Type", "RSD_before", "RSD_after", "delta_RSD"))
  expect_equal(nrow(comparison), 4L)
})

testthat::test_that("export_outliers_xlsx reuses the supplied detector result", {
  testthat::skip_if_not_installed("openxlsx")
  
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  hotelling_res <- make_cached_hotelling_result()
  
  testthat::local_mocked_bindings(
    detect_outliers = function(...) {
      stop("The supplied detector result should be reused.")
    }
  )
  
  session <- shiny::MockShinySession$new()
  on.exit(session$close(), add = TRUE)
  
  shiny::withReactiveDomain(
    session,
    export_outliers_xlsx(
      p = list(),
      d = list(filtered_corrected = list(df_no_mv = data.frame())),
      file = file,
      hotelling_res = hotelling_res
    )
  )
  
  testthat::expect_identical(
    openxlsx::getSheetNames(file),
    c("Potential Outlier Samples", "Study Extreme Values", "QC Extreme Values")
  )
  
  samples <- openxlsx::read.xlsx(
    file,
    sheet = "Potential Outlier Samples",
    startRow = 3
  )
  testthat::expect_named(
    samples,
    c("Sample", "Class", "Squared_Mahalanobis_distance")
  )
  testthat::expect_equal(samples$Sample, "S1")
  testthat::expect_equal(samples$Class, "Case")
  testthat::expect_equal(samples$Squared_Mahalanobis_distance, 12)
  
  expected_value_columns <- c(
    "Sample",
    "Class",
    "Metabolite",
    "Method",
    "Test_statistic",
    "Adjusted_p_value",
    "P_value_cutoff",
    "Squared_Mahalanobis_distance"
  )
  
  study <- openxlsx::read.xlsx(
    file,
    sheet = "Study Extreme Values",
    startRow = 3
  )
  qc <- openxlsx::read.xlsx(
    file,
    sheet = "QC Extreme Values",
    startRow = 3
  )
  
  testthat::expect_named(study, expected_value_columns)
  testthat::expect_named(qc, expected_value_columns)
  
  testthat::expect_equal(study$Sample, "S1")
  testthat::expect_equal(study$Class, "Case")
  testthat::expect_equal(study$Metabolite, "met_a")
  testthat::expect_equal(study$Method, "Grubbs-type")
  testthat::expect_equal(study$Test_statistic, 6.1234)
  testthat::expect_equal(study$Adjusted_p_value, 0.002)
  testthat::expect_equal(study$P_value_cutoff, 0.005)
  
  testthat::expect_equal(qc$Sample, "S2")
  testthat::expect_equal(qc$Class, "QC")
  testthat::expect_equal(qc$Metabolite, "met_a")
  testthat::expect_equal(qc$Method, "ROUT")
  testthat::expect_equal(qc$Test_statistic, -9.5678)
  testthat::expect_equal(qc$Adjusted_p_value, 0.0005)
  testthat::expect_equal(qc$P_value_cutoff, 0.001)
  
  # A QC candidate is exported even when its sample is inside the ellipse.
  testthat::expect_equal(qc$Squared_Mahalanobis_distance, 1)
  
  for (sheet in c("Study Extreme Values", "QC Extreme Values")) {
    description <- openxlsx::read.xlsx(
      file,
      sheet = sheet,
      rows = 1,
      colNames = FALSE
    )
    description_text <- paste(unlist(description, use.names = FALSE), collapse = " ")
    
    testthat::expect_match(description_text, "fewer than 8", fixed = TRUE)
    testthat::expect_match(description_text, "Grubbs-type", fixed = TRUE)
    testthat::expect_match(description_text, "Bonferroni-adjusted", fixed = TRUE)
    testthat::expect_match(description_text, "ROUT (Q = 0.1%)", fixed = TRUE)
    testthat::expect_match(description_text, "Flagged values were retained", fixed = TRUE)
  }
})

test_that("export_corr_xlsx preserves correlation workbook labels", {
  all_corr <- list(
    raw = data.frame(col1 = "met_a", col2 = "met_b", cor = 0.9, n_complete = 4),
    corrected = data.frame(col1 = "met_a", col2 = "met_b", cor = 0.8, n_complete = 4),
    transformed = data.frame(col1 = "met_a", col2 = "met_b", cor = 0.7, n_complete = 4),
    transformed_included = TRUE
  )
  file <- tempfile(fileext = ".xlsx")
  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(session, export_corr_xlsx(all_corr, file = file))

  expect_setequal(openxlsx::getSheetNames(file), c("Raw Data", "Corrected Data", "Transformed Data"))

  transformed <- openxlsx::read.xlsx(
    file,
    sheet = "Transformed Data",
    startRow = 3,
    sep.names = " "
  )
  expect_equal(names(transformed), c("Metabolite 1", "Metabolite 2", "Pearson's r", "n_complete"))
  expect_equal(transformed[["Pearson's r"]], 0.7)
})
