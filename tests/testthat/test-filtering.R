# test-filtering.R

# helper to build a small dataset
make_df <- function() {
  data.frame(
    sample = paste0("s", 1:6),
    batch = 1L,
    class = c("QC", "sample", "sample", "QC", "sample", "QC"),
    order = 1:6,
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
    batch = 1L,
    class = c("QC", "QC", "QC", "sample", "sample", "sample"),
    order = 1:6,
    A = c(qc_vals(tgt[1]), 10, 20, 30),
    B = c(qc_vals(tgt[2]), 5, 5, 5),
    C = c(qc_vals(tgt[3]), 1, 2, 3),
    D = c(qc_vals(tgt[4]), 7, 8, 9),
    check.names = FALSE
  )
}

test_that("filter_by_missing applies separate study and QC filtering", {
  df <- make_df()
  metab_cols <- c("A", "B", "C", "D")
  
  out <- filter_by_missing(
    df,
    metab_cols,
    mv_cutoff = 50
  )
  
  expect_named(
    out,
    c(
      "df",
      "mv_cutoff",
      "qc_mv_cutoff",
      "filter_rule",
      "mv_removed_cols",
      "study_mv_removed_cols",
      "qc_mv_removed_cols",
      "qc_missing_mets",
      "class_metab_all_missing"
    )
  )
  
  expect_equal(out$mv_cutoff, 50)
  expect_equal(out$qc_mv_cutoff, Inf)
  expect_identical(out$filter_rule, "any")
  
  # Study-sample missing-like percentages:
  # A: sample = 66.67% -> removed by study threshold
  # B: sample = 66.67% -> removed by study threshold
  # C: sample = 0%     -> retained
  # D: sample = 100%   -> removed by study threshold
  #
  # The default QC cutoff is Inf, so QC missingness does not remove
  # any additional metabolites.
  expect_setequal(
    names(out$df),
    c("sample", "batch", "class", "order", "C")
  )
  
  expect_setequal(
    out$study_mv_removed_cols,
    c("A", "B", "D")
  )
  
  expect_identical(
    out$qc_mv_removed_cols,
    character(0)
  )
  
  expect_setequal(
    out$mv_removed_cols,
    c("A", "B", "D")
  )
  
  # Among retained metabolites, C has one missing-like QC value because
  # its QC values are 0, 3, and 5.
  expect_identical(out$qc_missing_mets, "C")
  
  expect_s3_class(
    out$class_metab_all_missing,
    "data.frame"
  )
  
  expect_named(
    out$class_metab_all_missing,
    c("class", "metabolite", "n_rows_in_class")
  )
  
  expect_equal(
    nrow(out$class_metab_all_missing),
    0L
  )
})

test_that("filter_by_missing combines study and QC threshold removals", {
  df <- make_df()
  metab_cols <- c("A", "B", "C", "D")
  
  out <- filter_by_missing(
    df,
    metab_cols,
    mv_cutoff = 33.33,
    qc_mv_cutoff = 33.33
  )
  
  # Study filtering:
  # A: sample = 66.67% -> removed
  # B: sample = 66.67% -> removed
  # C: sample = 0%     -> retained by study filtering
  # D: sample = 100%   -> removed
  expect_setequal(
    out$study_mv_removed_cols,
    c("A", "B", "D")
  )
  
  # QC filtering:
  # A: QC = 0%         -> retained by QC filtering
  # B: QC = 33.333...% -> removed
  # C: QC = 33.333...% -> removed
  # D: QC = 100%       -> removed
  expect_setequal(
    out$qc_mv_removed_cols,
    c("B", "C", "D")
  )
  
  # The union of the study and QC removal lists contains all metabolites.
  expect_setequal(
    out$mv_removed_cols,
    c("A", "B", "C", "D")
  )
  
  expect_setequal(
    names(out$df),
    c("sample", "batch", "class", "order")
  )
  
  expect_identical(
    out$qc_missing_mets,
    character(0)
  )
  
  expect_equal(
    nrow(out$class_metab_all_missing),
    0L
  )
})

test_that("filter_by_missing reports retained class-metabolite pairs with all missing-like values", {
  df <- data.frame(
    sample = paste0("s", 1:6),
    batch = 1L,
    class = c(
      "QC",
      "QC",
      "sample",
      "sample",
      "sample",
      "sample"
    ),
    order = 1:6,
    A = c(1, 2, NA, NA, NA, NA),
    B = c(1, 2, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )
  
  out <- filter_by_missing(
    df,
    metab_cols = c("A", "B"),
    mv_cutoff = 100
  )
  
  expect_setequal(
    names(out$df),
    c("sample", "batch", "class", "order", "A", "B")
  )
  
  expect_identical(
    out$study_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$qc_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$mv_removed_cols,
    character(0)
  )
  
  expect_equal(
    nrow(out$class_metab_all_missing),
    1L
  )
  
  expect_identical(
    out$class_metab_all_missing$class,
    "sample"
  )
  
  expect_identical(
    out$class_metab_all_missing$metabolite,
    "A"
  )
  
  expect_identical(
    out$class_metab_all_missing$n_rows_in_class,
    4L
  )
})

test_that("filter_by_missing supports any and all study-class rules", {
  df <- data.frame(
    sample = paste0("s", 1:6),
    batch = 1L,
    class = c(
      "QC",
      "group1",
      "group1",
      "group2",
      "group2",
      "QC"
    ),
    order = 1:6,
    A = c(1, NA, NA, 2, 3, 1),
    B = c(1, NA, NA, NA, NA, 1),
    stringsAsFactors = FALSE
  )
  
  out_any <- filter_by_missing(
    df,
    metab_cols = c("A", "B"),
    mv_cutoff = 50,
    qc_mv_cutoff = Inf,
    filter_rule = "any"
  )
  
  out_all <- filter_by_missing(
    df,
    metab_cols = c("A", "B"),
    mv_cutoff = 50,
    qc_mv_cutoff = Inf,
    filter_rule = "all"
  )
  
  # A exceeds the threshold only in group1.
  # B exceeds the threshold in both study groups.
  expect_setequal(
    out_any$study_mv_removed_cols,
    c("A", "B")
  )
  
  expect_identical(
    out_all$study_mv_removed_cols,
    "B"
  )
  
  expect_setequal(
    names(out_any$df),
    c("sample", "batch", "class", "order")
  )
  
  expect_setequal(
    names(out_all$df),
    c("sample", "batch", "class", "order", "A")
  )
  
  expect_identical(out_any$filter_rule, "any")
  expect_identical(out_all$filter_rule, "all")
})

test_that("filter_by_missing errors if class column is absent", {
  df <- data.frame(
    sample = paste0("s", 1:3),
    order = 1:3,
    A = c(1, NA, 3),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    filter_by_missing(
      df,
      metab_cols = "A",
      mv_cutoff = 50
    ),
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
  
  out <- filter_by_missing(
    df,
    metab_cols = character(0),
    mv_cutoff = 50
  )
  
  expect_named(
    out,
    c(
      "df",
      "mv_cutoff",
      "qc_mv_cutoff",
      "filter_rule",
      "mv_removed_cols",
      "study_mv_removed_cols",
      "qc_mv_removed_cols",
      "qc_missing_mets",
      "class_metab_all_missing"
    )
  )
  
  expect_equal(out$df, df)
  expect_equal(out$mv_cutoff, 50)
  expect_equal(out$qc_mv_cutoff, Inf)
  expect_identical(out$filter_rule, "any")
  
  expect_identical(
    out$mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$study_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$qc_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$qc_missing_mets,
    character(0)
  )
  
  expect_equal(
    nrow(out$class_metab_all_missing),
    0L
  )
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
  
  out <- filter_by_missing(
    df,
    metab_cols = "A",
    mv_cutoff = 50
  )
  
  expect_setequal(
    names(out$df),
    c("sample", "batch", "class", "order", "A")
  )
  
  expect_identical(
    out$study_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$qc_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$qc_missing_mets,
    "A"
  )
})

test_that("filter_by_missing retains metabolites while study cutoff is NULL", {
  df <- data.frame(
    sample = paste0("s", 1:3),
    batch = 1L,
    class = c("QC", "sample", "QC"),
    order = 1:3,
    A = c(1, NA, 3),
    stringsAsFactors = FALSE
  )
  
  out <- filter_by_missing(
    df,
    metab_cols = "A",
    mv_cutoff = NULL
  )
  
  expect_setequal(
    names(out$df),
    c("sample", "batch", "class", "order", "A")
  )
  
  expect_null(out$mv_cutoff)
  expect_equal(out$qc_mv_cutoff, Inf)
  
  expect_identical(
    out$study_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$qc_mv_removed_cols,
    character(0)
  )
  
  expect_identical(
    out$mv_removed_cols,
    character(0)
  )
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

test_that("percent distance filtering returns flagged and retained metabolites", {
  df <- data.frame(
    sample = paste0("s", 1:6),
    batch = 1L,
    class = c("QC", "QC", "QC", "sample", "sample", "sample"),
    order = 1:6,
    near_qc = c(10, 10, 10, 11, 10, 9),
    far_qc = c(10, 10, 10, 30, 30, 30),
    zero_qc = c(0, 0, 0, 1, 1, 1),
    stringsAsFactors = FALSE
  )

  stats <- get_metabs_pct_diff_vs_qc_average(
    df,
    percent_threshold = 100,
    return_stats = TRUE
  )
  removed <- remove_metabs_pct_diff_vs_qc_average(
    df,
    percent_threshold = 100,
    return_result = TRUE
  )

  expect_named(stats, c(
    "metabolite",
    "sample_mean",
    "qc_mean",
    "percent_distance_from_qc_average",
    "flagged"
  ))
  expect_setequal(stats$metabolite[stats$flagged], c("far_qc", "zero_qc"))
  expect_setequal(removed$removed_metabolites, c("far_qc", "zero_qc"))
  expect_setequal(names(removed$df), c("sample", "batch", "class", "order", "near_qc"))
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
