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
    used_in_fit = c(TRUE, FALSE),
    check.names = FALSE
  )

  list(
    data = data,
    extreme_values = data.frame(
      sample = "S1",
      class = "Case",
      metabolite = "met_a",
      z_global = 4,
      abs_z_global = 4,
      z_class = 3.5,
      abs_z_class = 3.5,
      T2 = 12,
      check.names = FALSE
    ),
    pc_loadings = data.frame(metabolite = "met_a", PC1 = 0.8, PC2 = 0.2),
    z_global = data.frame(
      sample = c("S1", "S2"),
      batch = "B1",
      class = c("Case", "QC"),
      order = 1:2,
      met_a = c(4.1234, -0.1),
      check.names = FALSE
    ),
    z_class = data.frame(
      sample = c("S1", "S2"),
      batch = "B1",
      class = c("Case", "QC"),
      order = 1:2,
      met_a = c(3.5678, NA),
      check.names = FALSE
    ),
    pca_plot = ggplot2::ggplot()
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

test_that("export_outliers_xlsx reuses supplied Hotelling result", {
  file <- tempfile(fileext = ".xlsx")
  hotelling_res <- make_cached_hotelling_result()

  testthat::local_mocked_bindings(
    detect_hotelling_nonqc_dual_z = function(...) {
      stop("Hotelling result should be reused")
    }
  )

  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(
    session,
    export_outliers_xlsx(
      p = list(),
      d = list(filtered_corrected = list(df_no_mv = data.frame())),
      file = file,
      hotelling_res = hotelling_res
    )
  )

  samples <- openxlsx::read.xlsx(file, sheet = "Samples Outside Ellipse", startRow = 3)
  extremes <- openxlsx::read.xlsx(file, sheet = "Potential Extreme Values", startRow = 3)

  expect_equal(samples$Sample, "S1")
  expect_equal(extremes$Sample, "S1")
  expect_equal(extremes$Metabolite, "met_a")
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
