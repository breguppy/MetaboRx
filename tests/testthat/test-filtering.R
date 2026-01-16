# test-filtering.R

# helper to build a small dataset
make_df <- function() {
  data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC", "sample", "sample", "QC", "sample", "QC"),
    order  = 1:6,
    A = c(1, NA, 3, 4, NA, 6),
    # 33.33% missing
    B = c(NA, NA, NA, 2, 3, 4),
    # 50% missing
    C = c(0, 1, 2, 3, 4, 5),
    # 0% missing
    D = c(NA, NA, NA, NA, NA, NA),
    # 100% missing
    stringsAsFactors = FALSE
  )
}
# toy data frame that yields desired QC RSDs for A,B,C,D
make_df_for_rsd <- function(tgt = c(15, 25, NA_real_, 60)) {
  # helper to craft QC rows with target RSD %
  qc_vals <- function(target_pct, m = 100) {
    s <- m * target_pct / 100
    c(m - s, m, m + s)
  }
  data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC", "QC", "QC", "sample", "sample", "sample"),
    order  = 1:6,
    A = c(qc_vals(tgt[1]), 10, 20, 30),
    B = c(qc_vals(tgt[2]), 5, 5, 5),
    C = c(qc_vals(tgt[3]), 1, 2, 3),
    D = c(qc_vals(tgt[4]), 7, 8, 9),
    check.names = FALSE
  )
}

test_that("filter_by_missing keeps <= cutoff and reports removed and QC-missing",
          {
            df <- make_df()
            metab_cols <- c("A", "B", "C", "D")
            
            out <- filter_by_missing(df, metab_cols, mv_cutoff = 50)
            
            # structure
            expect_named(out,
                         c("df", "mv_cutoff", "mv_removed_cols", "qc_missing_mets", "class_metab_all_missing"))
            expect_equal(out$mv_cutoff, 50)
            
            # A (33.33) and C (0) kept. B (50) kept since <=. D removed.
            expect_setequal(names(out$df),
                            c("sample", "batch", "class", "order", "A", "B", "C"))
            expect_setequal(out$mv_removed_cols, "D")
            
            # QC rows are orders 1,4,6. Check QC-missing among kept cols.
            # In QC rows: A has no NA, B has 1 missing QC value, C has none missing.
            expect_length(out$qc_missing_mets, 1)
          })

test_that("filter_by_missing with strict cutoff removes equals when using < cutoff",
          {
            df <- make_df()
            metab_cols <- c("A", "B", "C", "D")
            
            # simulate strict policy by post-filtering test: create expected sets
            out <- filter_by_missing(df, metab_cols, mv_cutoff = 33.33)
            # Only C (0%) kept since A is ~33.33% > cutoff numerically (due to float).
            expect_setequal(names(out$df), c("sample", "batch", "class", "order", "C"))
            expect_setequal(out$mv_removed_cols, c("A", "B", "D"))
          })

test_that("remove_imputed_from_corrected masks positions where raw is NA",
          {
            raw <- data.frame(x = c(1, NA, 3), y = c(NA, 2, 3))
            cor <- data.frame(x = c(10, 20, 30), y = c(40, 50, 60))
            
            out <- remove_imputed_from_corrected(raw, cor)
            
            expect_equal(out$x, c(10, NA, 30))
            expect_equal(out$y, c(NA, 50, 60))
            # original objects unchanged
            expect_equal(cor$y[1], 40)
          })

test_that("remove_imputed_from_corrected errors on shape mismatch", {
  raw <- data.frame(x = 1:3, y = 1:3)
  cor <- data.frame(x = 1:3)
  expect_error(remove_imputed_from_corrected(raw, cor), "same dimensions")
})

test_that("filter_by_qc_rsd keeps <= cutoff, removes NA and > cutoff", {
  df <- make_df_for_rsd()
  
  out <- filter_by_qc_rsd(
    df,
    df,
    rsd_cutoff = 25,
    remove_imputed = TRUE,
    metadata_cols = c("sample", "batch", "class", "order")
  )
  
  # keep A (15) and C (NA). remove B (25? kept since <= cutoff) and D (60).
  expect_setequal(names(out$df_mv),
                  c("sample", "batch", "class", "order", "A", "B"))
  expect_setequal(out$removed_metabolites_mv, c("C", "D"))
  expect_equal(out$rsd_cutoff, 25)
})

test_that("filter_by_qc_rsd can remove all metabolites", {
  df <- make_df_for_rsd(c(70, 80, 90, 100))
  
  out <- filter_by_qc_rsd(df, df, rsd_cutoff = 60, remove_imputed = TRUE)
  
  expect_identical(setdiff(names(out$df_mv), c("sample", "batch", "class", "order")), character(0))
  expect_setequal(out$removed_metabolites_mv, c("A", "B", "C", "D"))
})

test_that("metabolite_rsd returns targeted QC RSDs", {
  df  <- make_df_for_rsd()
  rsd <- metabolite_rsd(df)
  got <- setNames(rsd$RSD_QC, rsd$Metabolite)
  
  expect_equal(unname(got["A"]), 15, tolerance = 1e-12)
  expect_equal(unname(got["B"]), 25, tolerance = 1e-12)
  expect_true(is.na(got["C"]))
  expect_equal(unname(got["D"]), 60, tolerance = 1e-12)
})