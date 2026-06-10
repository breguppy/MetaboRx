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
