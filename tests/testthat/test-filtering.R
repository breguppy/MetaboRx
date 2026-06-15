# test-filtering.R

# helper to build a small dataset
make_df <- function() {
  data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC", "sample", "sample", "QC", "sample", "QC"),
    order  = 1:6,
    A = c(1, NA, 3, 4, NA, 6),
    # sample class has 66.67% missing-like
    B = c(NA, NA, NA, 2, 3, 4),
    # QC = 33.33% missing-like, sample = 66.67% missing-like
    C = c(0, 1, 2, 3, 4, 5),
    # QC = 33.33% missing-like because 0 is treated as missing-like
    D = c(NA, NA, NA, NA, NA, NA),
    # 100% missing-like
    stringsAsFactors = FALSE
  )
}

# toy data frame that yields desired QC RSDs for A,B,C,D
make_df_for_rsd <- function(tgt = c(15, 25, NA_real_, 60)) {
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

test_that("filter_by_missing applies class-wise missing filtering and reports outputs correctly", {
  df <- make_df()
  metab_cols <- c("A", "B", "C", "D")
  
  out <- filter_by_missing(df, metab_cols, mv_cutoff = 50)
  
  expect_named(
    out,
    c("df", "mv_cutoff", "mv_removed_cols", "qc_missing_mets", "class_metab_all_missing")
  )
  expect_equal(out$mv_cutoff, 50)
  
  # Class-wise missing-like percentages:
  # A: QC = 0%, sample = 66.67% -> removed
  # B: QC = 33.33%, sample = 66.67% -> removed
  # C: QC = 33.33%, sample = 0% -> kept
  # D: 100% -> removed
  expect_setequal(names(out$df), c("sample", "batch", "class", "order", "C"))
  expect_setequal(out$mv_removed_cols, c("A", "B", "D"))
  
  # Among retained metabolites, QC rows for C are 0, 3, 5, so C has a QC missing-like value
  expect_identical(out$qc_missing_mets, "C")
  
  # No retained class-metabolite pair is entirely missing-like
  expect_s3_class(out$class_metab_all_missing, "data.frame")
  expect_named(out$class_metab_all_missing, c("class", "metabolite", "n_rows_in_class"))
  expect_equal(nrow(out$class_metab_all_missing), 0L)
})

test_that("filter_by_missing removes metabolites when any class exceeds the cutoff", {
  df <- make_df()
  metab_cols <- c("A", "B", "C", "D")
  
  out <- filter_by_missing(df, metab_cols, mv_cutoff = 33.33)
  
  # A: sample = 66.67% -> removed
  # B: sample = 66.67% -> removed
  # C: QC = 33.333...% which is > 33.33 -> removed
  # D: removed
  expect_setequal(names(out$df), c("sample", "batch", "class", "order"))
  expect_setequal(out$mv_removed_cols, c("A", "B", "C", "D"))
  expect_identical(out$qc_missing_mets, character(0))
  expect_equal(nrow(out$class_metab_all_missing), 0L)
})

test_that("filter_by_missing reports retained class-metabolite pairs with all missing-like values", {
  df <- data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC", "QC", "sample", "sample", "sample", "sample"),
    order  = 1:6,
    A = c(1, 2, NA, NA, NA, NA),   # retained at cutoff 100, but all missing-like in sample class
    B = c(1, 2, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )
  
  out <- filter_by_missing(df, metab_cols = c("A", "B"), mv_cutoff = 100)
  
  expect_setequal(names(out$df), c("sample", "batch", "class", "order", "A", "B"))
  expect_equal(nrow(out$class_metab_all_missing), 1L)
  expect_identical(out$class_metab_all_missing$class, "sample")
  expect_identical(out$class_metab_all_missing$metabolite, "A")
  expect_identical(out$class_metab_all_missing$n_rows_in_class, 4L)
})

test_that("filter_by_missing errors if class column is absent", {
  df <- data.frame(
    sample = paste0("s", 1:3),
    order = 1:3,
    A = c(1, NA, 3),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    filter_by_missing(df, metab_cols = "A", mv_cutoff = 50),
    "`df` must contain a 'class' column.",
    fixed = TRUE
  )
})

test_that("filter_by_missing handles empty metabolite sets", {
  df <- data.frame(
    sample = paste0("s", 1:3),
    batch = 1L,
    class = c("QC", "sample", "QC"),
    order = 1:3,
    stringsAsFactors = FALSE
  )

  out <- filter_by_missing(df, metab_cols = character(0), mv_cutoff = 50)

  expect_named(
    out,
    c("df", "mv_cutoff", "mv_removed_cols", "qc_missing_mets", "class_metab_all_missing")
  )
  expect_equal(out$df, df)
  expect_identical(out$mv_removed_cols, character(0))
  expect_identical(out$qc_missing_mets, character(0))
  expect_equal(nrow(out$class_metab_all_missing), 0L)
})

test_that("filter_by_missing keeps matrix shape with one metabolite and one class", {
  df <- data.frame(
    sample = paste0("s", 1:3),
    batch = 1L,
    class = c("QC", "QC", "QC"),
    order = 1:3,
    A = c(1, NA, 3),
    stringsAsFactors = FALSE
  )

  out <- filter_by_missing(df, metab_cols = "A", mv_cutoff = 50)

  expect_setequal(names(out$df), c("sample", "batch", "class", "order", "A"))
  expect_identical(out$mv_removed_cols, character(0))
  expect_identical(out$qc_missing_mets, "A")
})

test_that("filter_by_missing retains metabolites while cutoff is NULL", {
  df <- data.frame(
    sample = paste0("s", 1:3),
    batch = 1L,
    class = c("QC", "sample", "QC"),
    order = 1:3,
    A = c(1, NA, 3),
    stringsAsFactors = FALSE
  )

  out <- filter_by_missing(df, metab_cols = "A", mv_cutoff = NULL)

  expect_setequal(names(out$df), c("sample", "batch", "class", "order", "A"))
  expect_null(out$mv_cutoff)
  expect_identical(out$mv_removed_cols, character(0))
})

test_that("detect_blank_threshold returns expected vectorized threshold table", {
  df <- data.frame(
    sample = paste0("s", 1:4),
    batch = 1L,
    class = c("QC", "sample", "sample", "QC"),
    order = 1:4,
    met_high = c(10, 2, 2, 10),
    met_low = c(2, 10, 10, 2),
    ISTD_low = c(2, 10, 10, 2),
    stringsAsFactors = FALSE
  )
  blank_df <- data.frame(
    sample = paste0("b", 1:2),
    batch = 1L,
    class = "blank",
    order = 5:6,
    met_high = c(1, 1),
    met_low = c(1, 1),
    ISTD_low = c(1, 1),
    stringsAsFactors = FALSE
  )

  out <- detect_blank_threshold(
    df = df,
    blank_df = blank_df,
    metab_cols = c("met_high", "met_low", "ISTD_low"),
    threshold = 3
  )

  expect_named(out, c(
    "blank_means",
    "qc_means",
    "below_blank_threshold",
    "below_blank_threshold_ex_ISTD",
    "threshold_table"
  ))
  expect_setequal(out$below_blank_threshold, c("met_low", "ISTD_low"))
  expect_identical(out$below_blank_threshold_ex_ISTD, "met_low")
  expect_named(out$threshold_table, c(
    "metabolite",
    "blank_mean",
    "qc_mean",
    "threshold_value",
    "eligible",
    "below_blank_threshold",
    "internal_standard"
  ))
})

test_that("remove_imputed_from_corrected masks positions where raw is NA", {
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

test_that("filter_by_qc_rsd keeps <= cutoff and removes NA and > cutoff", {
  df <- make_df_for_rsd()
  
  out <- filter_by_qc_rsd(
    df,
    df,
    rsd_cutoff = 25,
    remove_imputed = TRUE,
    metadata_cols = c("sample", "batch", "class", "order")
  )
  
  # keep A (15) and B (25 <= cutoff); remove C (NA) and D (60)
  expect_setequal(
    names(out$df_mv),
    c("sample", "batch", "class", "order", "A", "B")
  )
  expect_setequal(out$removed_metabolites_mv, c("C", "D"))
  expect_equal(out$rsd_cutoff, 25)
})

test_that("filter_by_qc_rsd can remove all metabolites", {
  df <- make_df_for_rsd(c(70, 80, 90, 100))
  
  out <- filter_by_qc_rsd(df, df, rsd_cutoff = 60, remove_imputed = TRUE)
  
  expect_identical(
    setdiff(names(out$df_mv), c("sample", "batch", "class", "order")),
    character(0)
  )
  expect_setequal(out$removed_metabolites_mv, c("A", "B", "C", "D"))
})

test_that("metabolite_rsd returns targeted QC RSDs", {
  df <- make_df_for_rsd()
  rsd <- metabolite_rsd(df)
  got <- setNames(rsd$RSD_QC, rsd$Metabolite)
  
  expect_equal(unname(got["A"]), 15, tolerance = 1e-12)
  expect_equal(unname(got["B"]), 25, tolerance = 1e-12)
  expect_true(is.na(got["C"]))
  expect_equal(unname(got["D"]), 60, tolerance = 1e-12)
})
