test_that("export_mv_xlsx writes vectorized missing value summaries", {
  cleaned_df <- data.frame(
    sample = c("QC1", "S1", "S2", "QC2"),
    batch = c("B1", "B1", "B2", "B2"),
    class = c("QC", "Case", "Case", "QC"),
    order = 1:4,
    met_a = c(NA, 0, 3, 4),
    met_b = c(2, 5, -1, 6),
    check.names = FALSE
  )

  file <- tempfile(fileext = ".xlsx")
  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(
    session,
    export_mv_xlsx(list(), list(cleaned = list(df = cleaned_df)), file = file)
  )

  metabolite <- openxlsx::read.xlsx(file, sheet = "Metabolite", startRow = 3)
  sample <- openxlsx::read.xlsx(file, sheet = "Sample", startRow = 3)
  class <- openxlsx::read.xlsx(file, sheet = "Class", startRow = 3)
  batch <- openxlsx::read.xlsx(file, sheet = "Batch", startRow = 3)
  class_met <- openxlsx::read.xlsx(file, sheet = "Class-Met Missing", startRow = 3)

  expect_equal(metabolite$metabolite, c("met_a", "met_b"))
  expect_equal(metabolite$sample_missing_count, c(1, 1))
  expect_equal(metabolite$qc_missing_count, c(1, 0))

  expect_equal(sample$sample, c("QC1", "S1", "S2"))
  expect_equal(sample$missing_count, c(1, 1, 1))
  expect_equal(sample$missing_pct, c(50, 50, 50))

  expect_equal(class$class, c("Case", "QC"))
  expect_equal(class$missing_count, c(2, 1))
  expect_equal(class$missing_pct, c(50, 25))

  expect_equal(batch$batch, c("B1", "B2"))
  expect_equal(batch$missing_count, c(2, 1))
  expect_equal(batch$missing_pct, c(50, 25))

  expect_equal(class_met$class, c("Case", "Case", "QC", "QC"))
  expect_equal(class_met$metabolite, c("met_a", "met_b", "met_a", "met_b"))
  expect_equal(class_met$missing_count, c(1, 1, 1, 0))
})
