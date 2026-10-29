# Regression coverage for pre-imputation QC scaling references.
qc_reference_fixture <- function() {
  data.frame(
    sample = paste0("reference_s", 1:5), batch = "A",
    class = c("QC", "sample", "QC", "sample", "QC"),
    order = 1:5, M1 = c(10, 20, 100, 40, NA_real_),
    M2 = c(20, 40, 200, 80, NA_real_)
  )
}

reference_imputed <- function(raw, replacement) {
  out <- raw
  out$M1[is.na(out$M1)] <- replacement
  out$M2[is.na(out$M2)] <- 2 * replacement
  out
}

testthat::test_that("reference alignment follows sample IDs rather than row order", {
  raw <- qc_reference_fixture()
  shuffled <- raw[c(5, 3, 1, 4, 2), ]
  aligned <- .align_qc_reference(raw, shuffled, c("M1", "M2"))
  testthat::expect_identical(aligned$sample, raw$sample)
  testthat::expect_equal(aligned$M1, raw$M1)
  testthat::expect_identical(.align_qc_reference(raw, NULL, "M1"), raw)
})

testthat::test_that("reference alignment rejects ambiguous or incompatible inputs", {
  raw <- qc_reference_fixture()
  duplicate <- raw
  duplicate$sample[2] <- duplicate$sample[1]
  missing_id <- raw
  missing_id$sample[1] <- NA_character_
  wrong_class <- raw
  wrong_class$class[1] <- "sample"
  nonnumeric <- raw
  nonnumeric$M1 <- as.character(nonnumeric$M1)

  testthat::expect_error(.align_qc_reference(raw, duplicate, "M1"), "unique")
  testthat::expect_error(.align_qc_reference(duplicate, raw, "M1"), "unique")
  testthat::expect_error(.align_qc_reference(raw, missing_id, "M1"), "nonmissing")
  testthat::expect_error(.align_qc_reference(raw, raw[-1, ], "M1"), "missing correction samples")
  testthat::expect_error(.align_qc_reference(raw, wrong_class, "M1"), "classes")
  testthat::expect_error(.align_qc_reference(raw, nonnumeric, "M1"), "numeric")
  testthat::expect_error(.align_qc_reference(raw, raw["M1"], "M1"), "must contain")
  testthat::expect_error(.align_qc_reference(raw, numeric(0), "M1"), "must contain")
})

testthat::test_that("finalizer ignores externally imputed QCs for both targets", {
  raw <- qc_reference_fixture()
  for (replacement in c(10, 55, 1000)) {
    input <- reference_imputed(raw, replacement)
    # Arbitrary corrected values deliberately distinguish observed/all-QC medians.
    corrected <- input
    corrected$M1 <- c(2, 3, 8, 4, 100)
    corrected$M2 <- c(4, 6, 16, 8, 200)
    for (scaling in c("one", "original")) {
      out <- .finalize_qc_correction(
        corrected, input, c("M1", "M2"), c(1, 3, 5), scaling,
        reference_df = raw[c(5, 4, 3, 2, 1), ]
      )
      target <- if (scaling == "one") c(1, 1) else c(55, 110)
      testthat::expect_equal(stats::median(out$M1[c(1, 3)]), target[1])
      testthat::expect_equal(stats::median(out$M2[c(1, 3)]), target[2])
      testthat::expect_gt(out$M1[5], 0) # Imputed values remain in output.
    }
  }
})

testthat::test_that("raw target includes observed QCs whose predictions failed", {
  raw <- qc_reference_fixture()
  raw$M1[5] <- 1000
  corrected <- raw
  corrected$M1 <- c(2, 3, 8, 4, NA_real_)
  out <- .finalize_qc_correction(
    corrected, raw, "M1", c(1, 3, 5), "original",
    reference_df = raw, impute_missing = FALSE
  )
  # Raw target median(10, 100, 1000) = 100, not median(10, 100) = 55.
  testthat::expect_equal(stats::median(out$M1[c(1, 3)]), 100)
  testthat::expect_true(all(is.finite(out$M1)))
})

testthat::test_that(
  "no observed positive QCs triggers a warning and no invented target",
  {
    raw <- qc_reference_fixture()
    raw$M1[c(1, 3, 5)] <- c(0, NA_real_, Inf)
    
    input <- reference_imputed(raw, 10)
    input$M1[!is.finite(input$M1) | input$M1 <= 0] <- 10
    
    out <- NULL
    
    testthat::expect_warning(
      out <- .finalize_qc_correction(
        out = input,
        raw_df = input,
        metab_cols = "M1",
        qcid = c(1, 3, 5),
        met_scaling = "original",
        reference_df = raw
      ),
      regexp = "no valid positive observed QC pairs"
    )
    
    testthat::expect_s3_class(out, "data.frame")
    testthat::expect_equal(out$M1, input$M1)
  }
)

testthat::test_that("None uses observed QCs for one and preserves input for original", {
  raw <- qc_reference_fixture()
  input <- reference_imputed(raw, 10)
  one <- correct_data(input, c("M1", "M2"), "None", "one", reference_df = raw)
  original <- correct_data(input, c("M1", "M2"), "None", "original", reference_df = raw)
  testthat::expect_equal(one$df$M1, input$M1 / 55)
  testthat::expect_equal(one$df$M2, input$M2 / 110)
  testthat::expect_identical(original$df, input)
})

for (method_name in c("LC", "LL", "LOESS", "RF")) {
  local({
    method <- method_name
    testthat::test_that(paste(method, "passes the observed reference through orchestration"), {
      if (method == "RF") testthat::skip_if_not_installed("randomForest")
      raw <- qc_reference_fixture()
      # Six QC injections avoid the RF low-QC fallback warning. Only the
      # first two QCs were observed; all remaining QCs are externally imputed.
      extra <- raw[rep(5, 6), ]
      extra$sample <- paste0("reference_s", 6:11)
      extra$order <- 6:11
      extra$class <- ifelse(extra$order %% 2 == 1, "QC", "sample")
      raw <- rbind(raw, extra)
      rownames(raw) <- NULL
      for (replacement in c(10, 1000)) {
        input <- reference_imputed(raw, replacement)
        input <- input[c(11:6, 5, 3, 2, 4, 1), ]
        reference <- raw[c(3, 1, 5, 2, 4, 11:6), ]
        for (scaling in c("one", "original")) {
          out <- correct_data(
            input, c("M1", "M2"), method, scaling,
            reference_df = reference
          )$df
          observed_rows <- match(raw$sample[c(1, 3)], out$sample)
          target <- if (scaling == "one") c(1, 1) else c(55, 110)
          testthat::expect_equal(stats::median(out$M1[observed_rows]), target[1], tolerance = 1e-8)
          testthat::expect_equal(stats::median(out$M2[observed_rows]), target[2], tolerance = 1e-8)
          testthat::expect_true(all(is.finite(as.matrix(out[c("M1", "M2")]))))
        }
      }
    })
  })
}
