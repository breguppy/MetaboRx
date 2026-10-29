# test-rf.R

testthat::skip_if_not_installed("randomForest")

# --- helpers to make tiny datasets -----------------------------------------

mk_df_rf_single <- function(n = 11) {
  # sorted by order; QC at first and last (odd positions are QC)
  stopifnot(n >= 5)
  data.frame(
    sample = paste0("s", 1:n),
    batch  = "A",
    class  = ifelse((1:n) %% 2 == 1, "QC", "sample"),
    order  = 1:n,
    M1 = 100 + 2*(1:n),
    M2 =  50 + 1*(1:n),
    check.names = FALSE
  )
}

mk_df_rf_batches_warn <- function() {
  # two 6-row batches; QC at 1 and 6 in each; only 2 QCs per batch → warnings
  one <- function(tag) data.frame(
    sample = paste0(tag, "_s", 1:6),
    batch  = tag,
    class  = c("QC","sample","sample","sample","sample","QC"),
    order  = 1:6,
    M1 = 200 + 3*(1:6),
    M2 = 100 + 2*(1:6),
    check.names = FALSE
  )
  rbind(one("A"), one("B"))
}

mk_df_rf_batches_ok <- function(n = 11) {
  # two batches; QC at first and last; many QCs (odd positions)
  stopifnot(n >= 7)
  one <- function(tag) data.frame(
    sample = paste0(tag, "_s", 1:n),
    batch  = tag,
    class  = ifelse((1:n) %% 2 == 1, "QC", "sample"),
    order  = 1:n,
    M1 = 300 + 2*(1:n),
    M2 = 150 + 1.5*(1:n),
    check.names = FALSE
  )
  rbind(one("A"), one("B"))
}
# --- .safe_baseline_predict --------------------------------------------------

test_that(".safe_baseline_predict returns flat 1s with <2 finite QCs", {
  qcid <- c(2, 5)
  qc_y <- c(NA_real_, NA_real_)
  newx <- 1:6
  out <- .safe_baseline_predict(qcid, qc_y, newx)
  expect_equal(out, rep(1, length(newx)))
})

test_that(".safe_baseline_predict does linear interp and extrap", {
  qcid <- c(1, 3, 5)
  qc_y <- c(10, 20, 30)
  # interpolation at 2 and 4
  pred <- .safe_baseline_predict(qcid, qc_y, 1:5)
  expect_equal(pred[2], 15)
  expect_equal(pred[4], 25)
  # extrapolation (rule=2)
  expect_equal(.safe_baseline_predict(qcid, qc_y, 0), 10)
  expect_equal(.safe_baseline_predict(qcid, qc_y, 6), 30)
})

# --- .safe_rf_predict --------------------------------------------------------

test_that(".safe_rf_predict falls back to baseline when QC y not unique", {
  set.seed(1)
  qcid <- c(1, 3, 5, 7, 9)
  qcy  <- rep(100, length(qcid))  # all identical → baseline
  newx <- 1:10
  base <- .safe_baseline_predict(qcid, qcy, newx)
  got  <- .safe_rf_predict(qcid, qcy, newx, ntree = 50)
  expect_equal(got, base)
  expect_true(all(is.finite(got)))
})

test_that(".safe_rf_predict returns numeric finite vector with enough QCs", {
  set.seed(1)
  qcid <- c(1, 3, 5, 7, 9)
  qcy  <- 100 + 2*qcid + rnorm(length(qcid), 0, 0.1)
  newx <- 1:10
  got  <- .safe_rf_predict(qcid, qcy, newx, ntree = 100)
  expect_type(got, "double")
  expect_length(got, length(newx))
  expect_false(any(!is.finite(got)))
})

# --- rf_correction -----------------------------------------------------------

test_that("rf_correction returns finite, non-negative, no-NA; preserves shape", {
  df <- mk_df_rf_single(n = 11)
  out <- rf_correction(df, metab_cols = c("M1","M2"), ntree = 100, seed = 42, min_qc = 5)
  
  expect_setequal(names(out), names(df))
  mets <- as.matrix(out[c("M1","M2")])
  expect_false(any(is.na(mets)))
  expect_false(any(is.nan(mets)))
  expect_false(any(is.infinite(mets)))
  expect_true(all(mets >= 0))
  expect_equal(out$order, df$order)
})

test_that("rf_correction warns when < min_qc but still outputs valid data", {
  df <- mk_df_rf_single(n = 11)
  p <- testthat::evaluate_promise(
    rf_correction(df, metab_cols = c("M1","M2"), ntree = 50, seed = 1, min_qc = 12)
  )
  expect_true(any(grepl("Fewer than 12 QC samples; using baseline interpolation.", p$warnings)))
  out <- p$result
  mets <- as.matrix(out[c("M1","M2")])
  expect_false(any(!is.finite(mets)))
  expect_true(all(mets >= 0))
})

# --- bw_rf_correction --------------------------------------------------------

test_that("bw_rf_correction: batch with < min_qc uses baseline; output is clean", {
  df <- mk_df_rf_batches_warn()  # boundaries OK, min_qc triggers warning
  p <- testthat::evaluate_promise(
    bw_rf_correction(df, metab_cols = c("M1","M2"), ntree = 100, seed = 2, min_qc = 5)
  )
  # Both metabolites should warn for batch B
  expect_true(any(grepl("Skipping .* batch B - too few QC samples; using baseline", p$warnings)))
  out <- p$result
  mets <- as.matrix(out[c("M1","M2")])
  expect_false(any(!is.finite(mets)))
  expect_true(all(mets >= 0))
  # Check names and order preserved
  expect_setequal(names(out), names(df))
  expect_equal(out$order, df$order)
})

test_that("bw_rf_correction produces finite results when both batches have enough QCs", {
  df <- mk_df_rf_batches_ok(n = 11)
  out <- bw_rf_correction(df, metab_cols = c("M1","M2"), ntree = 100, seed = 3, min_qc = 5)
  mets <- as.matrix(out[c("M1","M2")])
  expect_false(any(!is.finite(mets)))
  expect_true(all(mets >= 0))
})

#--- QC boundary checks -----------------------------------------------------
  
  test_that("rf_correction errors if first/last are not QC", {
    df <- mk_df_rf_single(n = 10)              # last is sample
    expect_error(
      rf_correction(df, metab_cols = "M1"),
      "First and last samples must be QCs"
    )
    # fix boundary → should pass
    df_ok <- mk_df_rf_single(n = 11)           # last is QC
    expect_silent(
      rf_correction(df_ok, metab_cols = "M1", ntree = 50, seed = 1, min_qc = 5)
    )
  })

test_that("bw_rf_correction errors if a batch does not start/end with QC", {
  df <- mk_df_rf_batches_ok(n = 11)
  # break batch A end QC
  df$class[df$batch == "A" & df$order == max(df$order[df$batch == "A"])] <- "sample"
  expect_error(
    bw_rf_correction(df, metab_cols = "M1"),
    "Batch 'A' must start and end with QC"
  )
  # positive control: restore boundary → pass
  df_ok <- mk_df_rf_batches_ok(n = 11)
  expect_silent(
    bw_rf_correction(df_ok, metab_cols = "M1", ntree = 50, seed = 1, min_qc = 5)
  )
})