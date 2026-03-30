# test-rsd.R

# helpers
qc_triplet <- function(target_pct, m = 100) {
  if (is.na(target_pct)) return(rep(NA_real_, 3))
  s <- m * target_pct / 100
  c(m - s, m, m + s)
}
mk_df <- function() {
  data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC","QC","QC","sample","sample","sample"),
    order  = 1:6,
    A = c(qc_triplet(15), 10,20,30),
    B = c(qc_triplet(25),  5, 5, 5),
    C = c(qc_triplet(NA),  1, 2, 3),
    D = c(qc_triplet(60),  7, 8, 9),
    check.names = FALSE
  )
}

test_that("metabolite_rsd computes QC and NonQC RSDs and respects metadata casing", {
  df <- mk_df()
  names(df)[names(df) == "class"] <- "Class"  # case-insensitiv
  rsd <- metabolite_rsd(df, metadata_cols = c("sample","batch","class","order"))
  
  expect_setequal(names(rsd), c("Metabolite","RSD_QC","RSD_NonQC"))
  got_qc <- setNames(rsd$RSD_QC, rsd$Metabolite)
  expect_equal(unname(got_qc["A"]), 15, tolerance = 1e-12)
  expect_equal(unname(got_qc["B"]), 25, tolerance = 1e-12)
  expect_true(is.na(got_qc["C"]))
  expect_equal(unname(got_qc["D"]), 60, tolerance = 1e-12)
  
  # NonQC checks (quick sanity, not exact targets)
  expect_true(all(is.finite(rsd$RSD_NonQC[rsd$Metabolite %in% c("A","B","D")])))
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
  
  expect_setequal(names(out), c("class","Metabolite","Mean","SD","RSD"))
  # QC RSDs match targets
  qc <- subset(out, class == "QC")
  got <- setNames(qc$RSD, qc$Metabolite)
  expect_equal(unname(got["A"]), 15, tolerance = 1e-12)
  expect_equal(unname(got["B"]), 25, tolerance = 1e-12)
  expect_true(is.na(got["C"]))
  expect_equal(unname(got["D"]), 60, tolerance = 1e-12)
})

test_that("delta_rsd_stats works for wide schema (Metabolite, RSD_QC, RSD_NonQC)", {
  before <- data.frame(
    Metabolite = c("A","B","C"),
    RSD_QC = c(30, 20, 10),
    RSD_NonQC = c(40, 10, NA_real_)
  )
  after <- data.frame(
    Metabolite = c("A","B","C"),
    RSD_QC = c(10, 10, 20),
    RSD_NonQC = c(20,  5, NA_real_)
  )
  s <- delta_rsd_stats(before, after)
  expect_equal(s$avg_delta_qc,  mean(c(-20, -10, +10)))
  expect_equal(s$med_delta_qc,  stats::median(c(-20, -10, +10)))
  expect_equal(s$avg_delta_sample, mean(c(-20, -5), na.rm = TRUE))
  expect_equal(s$med_delta_sample, stats::median(c(-20, -5), na.rm = TRUE))
})

test_that("delta_rsd_stats works for long schema (class, Metabolite, RSD)", {
  # Before/After for QC and two sample classes
  before <- data.frame(
    class = c("QC","QC","QC", "A","A","A", "B","B","B"),
    Metabolite = rep(c("M1","M2","M3"), 3),
    RSD = c(30,20,10, 40,10,NA, 50, 5, 5)
  )
  after <- data.frame(
    class = c("QC","QC","QC", "A","A","A", "B","B","B"),
    Metabolite = rep(c("M1","M2","M3"), 3),
    RSD = c(10,10,20, 20, 5, NA, 30, 0, 5)
  )
  s <- delta_rsd_stats(before, after)
  # QC deltas: -20, -10, +10
  expect_equal(s$avg_delta_qc, mean(c(-20,-10,10)))
  expect_equal(s$med_delta_qc, stats::median(c(-20,-10,10)))
  # Sample deltas are now computed at the non-QC class-metabolite level:
  # A-M1 = -20, A-M2 = -5, A-M3 = NA
  # B-M1 = -20, B-M2 = -5, B-M3 = 0
  sample_deltas <- c(-20, -5, -20, -5, 0)
  
  expect_equal(s$avg_delta_sample, mean(sample_deltas))
  expect_equal(s$med_delta_sample, stats::median(sample_deltas))
})

test_that("delta_rsd_stats errors on unknown schema", {
  bad1 <- data.frame(Met = "A", RSD = 1)
  bad2 <- data.frame(Met = "A", RSD = 2)
  expect_error(delta_rsd_stats(bad1, bad2), "Unrecognized")
})
