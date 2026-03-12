# test-correcting.R

# --- helpers ---------------------------------------------------------------

mk_df_single <- function(n = 11) {
  stopifnot(n %% 2 == 1)  # ensures first/last row are QC
  data.frame(
    sample = paste0("s", seq_len(n)),
    batch  = "A",
    class  = ifelse(seq_len(n) %% 2 == 1, "QC", "sample"),
    order  = seq_len(n),
    M1     = 100 + 2 * seq_len(n),
    M2     =  50 + 1 * seq_len(n),
    check.names = FALSE
  )
}

mk_df_batches_ok <- function(n = 11) {
  stopifnot(n %% 2 == 1)
  
  one_batch <- function(tag) {
    data.frame(
      sample = paste0(tag, "_s", seq_len(n)),
      batch  = tag,
      class  = ifelse(seq_len(n) %% 2 == 1, "QC", "sample"),
      order  = seq_len(n),
      M1     = 300 + 2 * seq_len(n),
      M2     = 150 + 1.5 * seq_len(n),
      check.names = FALSE
    )
  }
  
  rbind(one_batch("A"), one_batch("B"))
}

meta_cols <- c("sample", "batch", "class", "order")
met_cols <- function(df) setdiff(names(df), meta_cols)

expect_clean_metabolites <- function(df, metab_cols) {
  mat <- as.matrix(df[metab_cols])
  testthat::expect_false(any(is.na(mat) | is.nan(mat) | is.infinite(mat)))
  testthat::expect_true(all(vapply(df[metab_cols], is.numeric, logical(1))))
  invisible(TRUE)
}

manual_seed_median <- function(dfs, metab_cols) {
  out <- dfs[[1]]
  out[metab_cols] <- Map(
    f = function(...) apply(cbind(...), 1, stats::median),
    dfs[[1]][metab_cols],
    dfs[[2]][metab_cols],
    dfs[[3]][metab_cols]
  )
  out
}

# --- Nadaraya-Watson / local constant --------------------------------------

testthat::test_that("correct_data LC returns clean, shaped output", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_single()
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "LC")
  
  testthat::expect_type(out, "list")
  testthat::expect_true(all(c("df", "str", "parameters") %in% names(out)))
  testthat::expect_equal(names(out$df), names(df))
  testthat::expect_equal(out$str, "local constant regression")
  testthat::expect_true(
    grepl("Nadaraya-Watson", out$parameters, fixed = TRUE) ||
      grepl("kernel-weighted mean", out$parameters, ignore.case = TRUE)
  )
  
  expect_clean_metabolites(out$df, met_cols(df))
})

# --- Local linear / local polynomial ---------------------------------------

testthat::test_that("correct_data LL returns clean, shaped output", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_single()
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "LL")
  
  testthat::expect_type(out, "list")
  testthat::expect_true(all(c("df", "str", "parameters") %in% names(out)))
  testthat::expect_equal(names(out$df), names(df))
  testthat::expect_equal(out$str, "local linear regression")
  testthat::expect_true(
    grepl("line created by nearby QC points", out$parameters, fixed = TRUE) ||
      grepl("local linear", out$parameters, ignore.case = TRUE)
  )
  
  expect_clean_metabolites(out$df, met_cols(df))
})

testthat::test_that("correct_data LOESS returns clean, shaped output", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_single()
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "LOESS")
  
  testthat::expect_type(out, "list")
  testthat::expect_true(all(c("df", "str", "parameters") %in% names(out)))
  testthat::expect_equal(names(out$df), names(df))
  testthat::expect_equal(out$str, "local polynomial regression")
  testthat::expect_true(
    grepl("degree 2", out$parameters, fixed = TRUE) ||
      grepl("local polynomials", out$parameters, ignore.case = TRUE)
  )
  
  expect_clean_metabolites(out$df, met_cols(df))
})

testthat::test_that("correct_data BW_LOESS returns clean, shaped output", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_batches_ok()
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "BW_LOESS")
  
  testthat::expect_type(out, "list")
  testthat::expect_true(all(c("df", "str", "parameters") %in% names(out)))
  testthat::expect_equal(names(out$df), names(df))
  testthat::expect_equal(out$str, "Batchwise LOESS")
  testthat::expect_true(
    grepl("degree 2", out$parameters, fixed = TRUE) ||
      grepl("local polynomials", out$parameters, ignore.case = TRUE)
  )
  
  expect_clean_metabolites(out$df, met_cols(df))
})

# --- Random forest ---------------------------------------------------------

testthat::test_that("correct_data RF equals median across the three RF seeds", {
  testthat::skip_if_not_installed("randomForest")
  
  df <- mk_df_single()
  seeds <- c(42, 31416, 272)
  
  dfs <- lapply(
    X = seeds,
    FUN = function(seed) {
      rf_correction(df, met_cols(df), ntree = 500, seed = seed)
    }
  )
  mdm <- manual_seed_median(dfs, met_cols(df))
  
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "RF")
  
  testthat::expect_equal(out$str, "Random Forest")
  testthat::expect_true(grepl("median value of the 3 models", out$parameters, fixed = TRUE))
  testthat::expect_equal(out$df[met_cols(df)], mdm[met_cols(df)], tolerance = 1e-8)
  testthat::expect_equal(names(out$df), names(df))
})

testthat::test_that("correct_data BW_RF equals median across batch-wise RF seeds", {
  testthat::skip_if_not_installed("randomForest")
  
  df <- mk_df_batches_ok()
  seeds <- c(42, 31416, 272)
  
  dfs <- lapply(
    X = seeds,
    FUN = function(seed) {
      bw_rf_correction(df, met_cols(df), ntree = 500, seed = seed)
    }
  )
  mdm <- manual_seed_median(dfs, met_cols(df))
  
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "BW_RF")
  
  ord_out <- order(out$df$batch, out$df$order, out$df$sample)
  ord_mdm <- order(mdm$batch, mdm$order, mdm$sample)
  
  testthat::expect_equal(out$str, "Batchwise Random Forest")
  testthat::expect_true(grepl("median value of the 3 models", out$parameters, fixed = TRUE))
  testthat::expect_equal(
    out$df[ord_out, met_cols(df), drop = FALSE],
    mdm[ord_mdm, met_cols(df), drop = FALSE],
    tolerance = 1e-8
  )
  testthat::expect_setequal(names(out$df), names(df))
})

# --- Direct local-constant implementation tests ----------------------------

testthat::test_that("nw_correction returns clean output with preserved metadata", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_single()
  out <- nw_correction(df, metab_cols = met_cols(df), span = 0.75)
  
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(names(out), names(df))
  testthat::expect_identical(out[meta_cols], df[meta_cols])
  
  expect_clean_metabolites(out, met_cols(df))
})

testthat::test_that("nw_correction preserves exact zeros in metabolite columns", {
  testthat::skip_if_not_installed("impute")
  
  df <- mk_df_single()
  df$M1[c(2, 6)] <- 0
  df$M2[c(4, 8)] <- 0
  
  out <- nw_correction(df, metab_cols = met_cols(df), span = 0.75)
  
  testthat::expect_identical(out$M1[c(2, 6)], c(0, 0))
  testthat::expect_identical(out$M2[c(4, 8)], c(0, 0))
  expect_clean_metabolites(out, met_cols(df))
})

testthat::test_that(".safe_nw_predict_x returns finite positive predictions on simple input", {
  qc_x <- c(1, 3, 5, 7, 9)
  qc_y <- c(10, 12, 14, 16, 18)
  newx <- 1:9
  
  pred <- .safe_nw_predict_x(qc_x = qc_x, qc_y = qc_y, newx = newx, span = 0.75)
  
  testthat::expect_true(is.numeric(pred))
  testthat::expect_equal(length(pred), length(newx))
  testthat::expect_false(any(is.na(pred) | is.nan(pred) | is.infinite(pred)))
  testthat::expect_true(all(pred > 0))
})

testthat::test_that(".safe_nw_predict_x falls back cleanly when too few QC points are available", {
  qc_x <- c(1)
  qc_y <- c(10)
  newx <- c(1, 2, 3, 4)
  
  pred <- .safe_nw_predict_x(qc_x = qc_x, qc_y = qc_y, newx = newx, span = 0.75)
  
  testthat::expect_equal(pred, rep(1, length(newx)))
})

# --- Generic properties ----------------------------------------------------

testthat::test_that("correct_data preserves metadata columns and types for single-batch methods", {
  testthat::skip_if_not_installed("impute")
  testthat::skip_if_not_installed("randomForest")
  
  df <- mk_df_single()
  methods <- c("LC", "LL", "LOESS", "RF")
  
  for (method in methods) {
    out <- correct_data(df, metab_cols = met_cols(df), corMethod = method)
    testthat::expect_true(all(meta_cols %in% names(out$df)))
    testthat::expect_identical(out$df[meta_cols], df[meta_cols])
  }
})

testthat::test_that("correct_data preserves metadata columns for batch-wise methods", {
  testthat::skip_if_not_installed("impute")
  testthat::skip_if_not_installed("randomForest")
  
  df <- mk_df_batches_ok()
  
  out_bw_loess <- correct_data(df, metab_cols = met_cols(df), corMethod = "BW_LOESS")
  out_bw_rf <- correct_data(df, metab_cols = met_cols(df), corMethod = "BW_RF")
  
  testthat::expect_identical(out_bw_loess$df[meta_cols], df[meta_cols])
  testthat::expect_identical(out_bw_rf$df[meta_cols], df[meta_cols])
})

testthat::test_that("correct_data returns numeric non-negative metabolite values across all methods", {
  testthat::skip_if_not_installed("impute")
  testthat::skip_if_not_installed("randomForest")
  
  single_df <- mk_df_single()
  batch_df <- mk_df_batches_ok()
  
  single_methods <- c("LC", "LL", "LOESS", "RF")
  batch_methods <- c("BW_LOESS", "BW_RF")
  
  for (method in single_methods) {
    out <- correct_data(single_df, metab_cols = met_cols(single_df), corMethod = method)
    expect_clean_metabolites(out$df, met_cols(single_df))
    testthat::expect_true(all(as.matrix(out$df[met_cols(single_df)]) >= 0))
  }
  
  for (method in batch_methods) {
    out <- correct_data(batch_df, metab_cols = met_cols(batch_df), corMethod = method)
    expect_clean_metabolites(out$df, met_cols(batch_df))
    testthat::expect_true(all(as.matrix(out$df[met_cols(batch_df)]) >= 0))
  }
})
