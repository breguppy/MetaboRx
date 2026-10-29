# test-rsd.R

# helpers
qc_triplet <- function(target_pct, m = 100) {
  if (is.na(target_pct)) {
    return(rep(NA_real_, 3))
  }
  s <- m * target_pct / 100
  c(m - s, m, m + s)
}
mk_df <- function() {
  data.frame(
    sample = paste0("s", 1:6),
    batch = 1L,
    class = c("QC", "QC", "QC", "sample", "sample", "sample"),
    order = 1:6,
    A = c(qc_triplet(15), 10, 20, 30),
    B = c(qc_triplet(25), 5, 5, 5),
    C = c(qc_triplet(NA), 1, 2, 3),
    D = c(qc_triplet(60), 7, 8, 9),
    check.names = FALSE
  )
}

test_that("metabolite_rsd computes QC and NonQC RSDs and respects metadata casing", {
  df <- mk_df()
  names(df)[names(df) == "class"] <- "Class" # case-insensitiv
  rsd <- metabolite_rsd(df, metadata_cols = c("sample", "batch", "class", "order"))

  expect_setequal(names(rsd), c("Metabolite", "RSD_QC", "RSD_NonQC"))
  got_qc <- setNames(rsd$RSD_QC, rsd$Metabolite)
  expect_equal(unname(got_qc["A"]), 15, tolerance = 1e-12)
  expect_equal(unname(got_qc["B"]), 25, tolerance = 1e-12)
  expect_true(is.na(got_qc["C"]))
  expect_equal(unname(got_qc["D"]), 60, tolerance = 1e-12)

  # NonQC checks (quick sanity, not exact targets)
  expect_true(all(is.finite(rsd$RSD_NonQC[rsd$Metabolite %in% c("A", "B", "D")])))
})

test_that("metabolite_rsd errors without class or without numeric metabolites", {
  df <- mk_df()
  df_no_class <- df[, setdiff(names(df), "class")]
  expect_error(metabolite_rsd(df_no_class), "Expected a 'class'")

  df_nonnum <- df
  df_nonnum$A <- as.character(df_nonnum$A)
  df_nonnum$B <- as.character(df_nonnum$B)
  df_nonnum$C <- as.character(df_nonnum$C)
  df_nonnum$D <- as.character(df_nonnum$D)
  expect_error(metabolite_rsd(df_nonnum), "No numeric metabolite")
})

test_that("class_metabolite_rsd returns per-class RSD with expected values", {
  df <- mk_df()
  out <- class_metabolite_rsd(df)

  expect_setequal(names(out), c("class", "Metabolite", "Mean", "SD", "RSD"))
  # QC RSDs match targets
  qc <- subset(out, class == "QC")
  got <- setNames(qc$RSD, qc$Metabolite)
  expect_equal(unname(got["A"]), 15, tolerance = 1e-12)
  expect_equal(unname(got["B"]), 25, tolerance = 1e-12)
  expect_true(is.na(got["C"]))
  expect_equal(unname(got["D"]), 60, tolerance = 1e-12)
})

# --- Tolerance-aware RSD differences ---------------------------------------

test_that(".rsd_change_delta returns signed changes for valid pairs", {
  expect_equal(
    .rsd_change_delta(
      before = c(10, 20, 30),
      after = c(12, 17, 30)
    ),
    c(2, -3, 0)
  )
})

test_that(".rsd_change_delta handles absolute tolerance boundaries", {
  # Binary-exact values avoid ambiguity at the equality boundary.
  before <- rep(0, 7)
  after <- c(0, 0.0625, -0.0625, 0.125, -0.125, 0.25, -0.25)
  
  expect_equal(
    .rsd_change_delta(
      before,
      after,
      atol = 0.125,
      rtol = 0
    ),
    c(0, 0, 0, 0, 0, 0.25, -0.25)
  )
})

test_that(".rsd_change_delta applies relative and combined tolerances", {
  # Relative tolerance uses max(abs(before), abs(after)).
  expect_equal(
    .rsd_change_delta(7, 8, atol = 0, rtol = 0.125),
    0
  )
  expect_equal(
    .rsd_change_delta(8, 7, atol = 0, rtol = 0.125),
    0
  )
  expect_equal(
    .rsd_change_delta(6, 8, atol = 0, rtol = 0.125),
    2
  )
  
  # Threshold = 1 + 0.125 * 8 = 2: equality means unchanged.
  expect_equal(
    .rsd_change_delta(6, 8, atol = 1, rtol = 0.125),
    0
  )
})

test_that(".rsd_change_delta preserves nonzero changes with zero tolerance", {
  expect_equal(
    .rsd_change_delta(
      c(0, 0),
      c(1e-12, -1e-12),
      atol = 0,
      rtol = 0
    ),
    c(1e-12, -1e-12),
    tolerance = 0
  )
})

test_that(".rsd_change_delta excludes nonfinite pairs", {
  expect_identical(
    .rsd_change_delta(
      before = c(10, NA_real_, NaN, Inf, -Inf, 10, 10, 10),
      after = c(12, 10, 10, 10, 10, NA_real_, NaN, Inf)
    ),
    c(2, rep(NA_real_, 7))
  )
  
  # Finite inputs can also produce an overflowing difference.
  expect_identical(
    .rsd_change_delta(-1e308, 1e308),
    NA_real_
  )
})

test_that(".rsd_change_delta accepts empty numeric inputs", {
  expect_identical(
    .rsd_change_delta(numeric(0), numeric(0)),
    numeric(0)
  )
})

# --- Summaries of paired RSD changes ---------------------------------------

test_that(".summarize_rsd_changes reports signed changes and percentages", {
  result <- .summarize_rsd_changes(
    rsd_before = c(10, 20, 30, 40),
    rsd_after = c(12, 17, 30, 44)
  )
  
  expect_equal(
    result,
    data.frame(
      med_delta = 1,
      pct_increase = 50,
      pct_decrease = 25,
      pct_unchanged = 25
    )
  )
})

test_that(".summarize_rsd_changes treats inside and boundary changes as zero", {
  result <- .summarize_rsd_changes(
    rsd_before = rep(0, 6),
    rsd_after = c(0.0625, -0.0625, 0.125, -0.125, 0.25, -0.25),
    atol = 0.125,
    rtol = 0
  )
  
  expect_equal(result$med_delta, 0)
  expect_equal(result$pct_increase, 100 / 6)
  expect_equal(result$pct_decrease, 100 / 6)
  expect_equal(result$pct_unchanged, 400 / 6)
})

test_that(".summarize_rsd_changes honors relative tolerance", {
  result <- .summarize_rsd_changes(
    rsd_before = c(7, 8, 6, 8),
    rsd_after = c(8, 7, 8, 6),
    atol = 0,
    rtol = 0.125
  )
  
  expect_equal(
    result,
    data.frame(
      med_delta = 0,
      pct_increase = 25,
      pct_decrease = 25,
      pct_unchanged = 50
    )
  )
})

test_that(".summarize_rsd_changes excludes invalid pairs from denominators", {
  result <- .summarize_rsd_changes(
    rsd_before = c(10, 20, 30, NA_real_, NaN, Inf, 40, -Inf),
    rsd_after = c(12, 18, 30, 10, 10, 10, Inf, 10)
  )
  
  expect_equal(result$med_delta, 0)
  expect_equal(result$pct_increase, 100 / 3)
  expect_equal(result$pct_decrease, 100 / 3)
  expect_equal(result$pct_unchanged, 100 / 3)
})

test_that(".summarize_rsd_changes handles empty and entirely invalid inputs", {
  expected <- data.frame(
    med_delta = NA_real_,
    pct_increase = NA_real_,
    pct_decrease = NA_real_,
    pct_unchanged = NA_real_
  )
  
  expect_identical(
    .summarize_rsd_changes(numeric(0), numeric(0)),
    expected
  )
  expect_identical(
    .summarize_rsd_changes(
      c(NA_real_, NaN, Inf, -Inf),
      c(1, 2, 3, 4)
    ),
    expected
  )
})

test_that(".summarize_rsd_changes handles a single valid pair", {
  expect_equal(
    .summarize_rsd_changes(10, 12),
    data.frame(
      med_delta = 2,
      pct_increase = 100,
      pct_decrease = 0,
      pct_unchanged = 0
    )
  )
})

# --- Argument validation for both new helpers ------------------------------

test_that("RSD change helpers reject invalid paired inputs", {
  helpers <- list(.rsd_change_delta, .summarize_rsd_changes)
  
  invalid_pairs <- list(
    list("10", 12),
    list(10, "12"),
    list(TRUE, 12),
    list(list(10), 12),
    list(NULL, numeric(0)),
    list(c(10, 20), 12),
    list(numeric(0), 12)
  )
  
  for (helper in helpers) {
    for (pair in invalid_pairs) {
      expect_error(helper(pair[[1]], pair[[2]]))
    }
  }
})

test_that("RSD change helpers reject invalid tolerances", {
  helpers <- list(.rsd_change_delta, .summarize_rsd_changes)
  
  invalid_tolerances <- list(
    -1,
    NA_real_,
    NaN,
    Inf,
    -Inf,
    numeric(0),
    c(0, 1),
    "invalid"
  )
  
  for (helper in helpers) {
    for (value in invalid_tolerances) {
      expect_error(helper(10, 12, atol = value))
      expect_error(helper(10, 12, rtol = value))
    }
    
    # Validate arguments even when there are no observations.
    expect_error(helper(numeric(0), numeric(0), atol = -1))
    expect_error(helper(numeric(0), numeric(0), rtol = -1))
  }
})

# --- Report integration ----------------------------------------------------

test_that("RSD labels distinguish numerical noise from real changes", {
  labels <- .rsd_change_label(
    before = c(20, 20, 20, 20, NA_real_, Inf),
    after = c(20 + 1e-12, 21, 19, 20, 20, 20)
  )
  
  expect_identical(
    labels,
    c(
      "No Change",
      "Increased",
      "Decreased",
      "No Change",
      NA_character_,
      NA_character_
    )
  )
})

test_that("RSD report summaries agree with tolerance-aware pair summaries", {
  comparisons <- data.frame(
    Type = c("QC", "QC", "QC", "QC", "Samples", "Samples"),
    before = c(20, 20, 20, NA_real_, 30, 30),
    after = c(20 + 1e-12, 22, 18, 20, 30 + 1e-12, 33)
  )
  
  # Deliberately provide unadjusted deltas to verify the summary applies
  # tolerance itself rather than trusting existing differences.
  comparisons$delta <- comparisons$after - comparisons$before
  
  result <- .summarize_rsd_comparison(comparisons)
  
  for (type in c("QC", "Samples")) {
    rows <- comparisons$Type == type
    expected <- .summarize_rsd_changes(
      comparisons$before[rows],
      comparisons$after[rows]
    )
    actual <- result[result$Type == type, , drop = FALSE]
    
    expect_equal(actual$med_delta, expected$med_delta)
    expect_equal(actual$pct_increase, expected$pct_increase)
    expect_equal(actual$pct_decrease, expected$pct_decrease)
    expect_equal(actual$pct_unchanged, expected$pct_unchanged)
  }
})
