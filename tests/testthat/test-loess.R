# test-loess.R

set.seed(1)

mk_df_single <- function() {
  data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("sample", "QC", "sample", "sample", "sample", "QC"),
    order  = c(2, 1, 3, 4, 5, 6),
    M1 = c(110, 100, 120, 130, 140, 150),
    M2 = c(55, 50, 52, 60, 58, 65),
    check.names = FALSE
  )
}

mk_df_batches <- function() {
  df <- data.frame(
    sample = paste0("b", rep(1:2, each = 6), "_s", 1:6),
    batch  = rep(c("A", "B"), each = 6),
    class  = rep(c("QC", "sample", "sample", "sample", "sample", "QC"), 2),
    order  = rep(1:6, 2),
    M1 = c(100, 105, 110, 115, 120, 125, 200, 195, 190, 185, 180, 175),
    M2 = c(50, 48, 51, 49, 47, 46, 100, 102, 104, 106, 108, 110),
    check.names = FALSE
  )
  df$class[df$batch == "B" & df$order %in% 2:5] <- "sample"
  df
}

test_that("loess_correction errors if first/last are not QC", {
  df <- mk_df_single()
  df$class[df$order == 1] <- "sample"
  
  expect_error(
    loess_correction(df, metab_cols = c("M1", "M2"), degree = 2),
    "First and last samples must be QCs"
  )
})

test_that("loess_correction sorts by order and preserves columns", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_single()
  out <- loess_correction(df, metab_cols = c("M1", "M2"), degree = 2, min_qc = 2)
  
  expect_equal(out$order, sort(df$order))
  expect_setequal(names(out), names(df))
})

test_that("loess_correction returns finite non-negative values and centers QC near 1", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_single()
  expect_silent({
    out <- loess_correction(df, metab_cols = c("M1", "M2"), degree = 2, min_qc = 2)
  })
  
  mets <- as.matrix(out[c("M1", "M2")])
  qc <- out[out$class == "QC", c("M1", "M2"), drop = FALSE]
  
  expect_false(any(is.na(mets)))
  expect_false(any(is.nan(mets)))
  expect_false(any(is.infinite(mets)))
  expect_true(all(mets >= 0))
  expect_true(all(vapply(out[c("M1", "M2")], is.numeric, TRUE)))
  
  expect_lt(abs(stats::median(qc$M1) - 1), 0.35)
  expect_lt(abs(stats::median(qc$M2) - 1), 0.35)
})

test_that("bw_loess_correction errors if a batch does not start and end with QC", {
  df <- mk_df_batches()
  df$class[df$batch == "A" & df$order == 6] <- "sample"
  
  expect_error(
    bw_loess_correction(df, metab_cols = c("M1", "M2"), degree = 2),
    "must start and end with QC"
  )
})

test_that("bw_loess_correction warns when batches have too few QCs", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_batches()
  p <- testthat::evaluate_promise(
    bw_loess_correction(df, metab_cols = c("M1", "M2"), min_qc = 5, degree = 1)
  )
  
  expect_true(any(grepl("Skipping batch 'A'", p$warnings)))
  expect_true(any(grepl("Skipping batch 'B'", p$warnings)))
  expect_setequal(names(p$result), names(df))
})

test_that("bw_loess_correction returns finite non-negative values and centers QC near 1 by batch", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_batches()
  expect_silent({
    out <- bw_loess_correction(df, metab_cols = c("M1", "M2"), min_qc = 2, degree = 1)
  })
  
  mets <- as.matrix(out[c("M1", "M2")])
  qc_A <- subset(out, batch == "A" & class == "QC", select = c("M1", "M2"))
  qc_B <- subset(out, batch == "B" & class == "QC", select = c("M1", "M2"))
  
  expect_false(any(is.na(mets)))
  expect_false(any(is.nan(mets)))
  expect_false(any(is.infinite(mets)))
  expect_true(all(mets >= 0))
  expect_true(all(vapply(out[c("M1", "M2")], is.numeric, TRUE)))
  
  expect_lt(abs(stats::median(qc_A$M1) - 1), 0.35)
  expect_lt(abs(stats::median(qc_A$M2) - 1), 0.35)
  expect_lt(abs(stats::median(qc_B$M1) - 1), 0.35)
  expect_lt(abs(stats::median(qc_B$M2) - 1), 0.35)
})

test_that("loess_correction cleanup falls back to zero when a metabolite is entirely zero", {
  testthat::skip_if_not_installed("impute")
  
  df <- data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC", "sample", "sample", "sample", "sample", "QC"),
    order  = 1:6,
    M3     = c(0, 0, 0, 0, 0, 0),
    check.names = FALSE
  )
  
  out <- loess_correction(df, metab_cols = "M3", degree = 2, min_qc = 2)
  expect_true(all(out$M3 == 0))
})

test_that("bw_loess_correction cleanup uses smallest positive or zero fallback", {
  testthat::skip_if_not_installed("impute")
  
  df <- data.frame(
    sample = paste0("b", rep(c("A", "B"), each = 6), "_s", 1:6),
    batch  = rep(c("A", "B"), each = 6),
    class  = rep(c("QC", "sample", "sample", "sample", "sample", "QC"), 2),
    order  = rep(1:6, 2),
    M4     = c(1, 2, NA, 3, 4, 5, 0, 0, 0, 0, 0, 0),
    check.names = FALSE
  )
  
  out <- bw_loess_correction(df, metab_cols = "M4", min_qc = 2, degree = 1)
  
  expect_false(any(is.na(out$M4)))
  expect_true(all(out$M4 >= 0))
  
  mp <- suppressWarnings(min(out$M4[out$M4 > 0]))
  bvals <- out$M4[out$batch == "B"]
  
  if (is.finite(mp)) {
    expect_true(all(bvals %in% c(0, mp)))
  } else {
    expect_true(all(bvals == 0))
  }
})