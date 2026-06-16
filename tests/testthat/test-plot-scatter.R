make_scatter_df <- function() {
  data.frame(
    sample = paste0("S", seq_len(8)),
    batch = rep(c("B1", "B2"), each = 4),
    class = c("QC", "Sample", "Sample", "QC", "QC", "Sample", "Sample", "QC"),
    order = seq_len(8),
    met_a = c(10, 12, 13, 11, 10.5, 12.5, 13.5, 11.5),
    met_b = c(20, 21, 22, 20.5, 20.2, 21.2, 22.2, 20.7),
    check.names = FALSE
  )
}

test_that("scatter context keeps only selected metabolite and plotting columns", {
  raw_df <- make_scatter_df()
  corrected_df <- dplyr::mutate(raw_df, met_a = met_a / mean(met_a))

  context <- .scatter_prepare_context(raw_df, corrected_df, "met_a")

  testthat::expect_setequal(
    names(context$df_all),
    c("sample", "batch", "class", "order", "met_a", "type", "panel")
  )
  testthat::expect_false("met_b" %in% names(context$df_all))
  testthat::expect_equal(nrow(context$df_all), 2L * nrow(raw_df))
  testthat::expect_equal(levels(context$df_all$panel), c("Raw", "Corrected"))
})

test_that("RF scatter accepts precomputed context without changing visible labels", {
  raw_df <- make_scatter_df()
  corrected_df <- dplyr::mutate(raw_df, met_a = met_a / mean(met_a))
  context <- .scatter_prepare_context(raw_df, corrected_df, "met_a")

  plot <- met_scatter_rf(raw_df, corrected_df, "met_a", scatter_context = context)
  built <- ggplot2::ggplot_build(plot)

  point_rows <- vapply(
    built$data,
    function(layer_data) {
      if (all(c("x", "y", "colour") %in% names(layer_data))) {
        nrow(layer_data)
      } else {
        0L
      }
    },
    integer(1)
  )

  testthat::expect_s3_class(plot, "ggplot")
  testthat::expect_equal(plot[["labels"]][["title"]], "met_a")
  testthat::expect_equal(plot[["labels"]][["x"]], "Injection Order")
  testthat::expect_equal(plot[["labels"]][["y"]], "Intensity")
  testthat::expect_equal(sum(point_rows), 2L * nrow(raw_df))
})

test_that("LOESS scatter accepts precomputed context without changing visible labels", {
  raw_df <- make_scatter_df()
  corrected_df <- dplyr::mutate(raw_df, met_a = met_a / mean(met_a))
  context <- .scatter_prepare_context(raw_df, corrected_df, "met_a")

  plot <- met_scatter_loess(
    raw_df,
    corrected_df,
    "local linear regression",
    "met_a",
    scatter_context = context
  )

  testthat::expect_s3_class(plot, "ggplot")
  testthat::expect_equal(plot[["labels"]][["title"]], "met_a")
  testthat::expect_equal(plot[["labels"]][["x"]], "Injection Order")
  testthat::expect_equal(plot[["labels"]][["y"]], "Intensity")
})
