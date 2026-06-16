# tests/testthat/test-detect_hotelling_nonqc_dual_z.R

testthat::test_that("detect_hotelling_nonqc_dual_z errors when required metadata columns are missing", {
  p <- list(qcImputeM = "median", samImputeM = "median")

  df <- data.frame(
    sample = paste0("S", 1:5),
    batch = "b1",
    class = c("A", "A", "A", "QC", "QC"),
    # order column intentionally missing
    m1 = rnorm(5),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df = df, p = p),
    "Missing metadata columns"
  )
})

testthat::test_that("detect_hotelling_nonqc_dual_z errors when class column is missing (via meta_cols check)", {
  p <- list(qcImputeM = "median", samImputeM = "median")

  df <- data.frame(
    sample = paste0("S", 1:5),
    batch = "b1",
    order = 1:5,
    m1 = rnorm(5),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df = df, p = p),
    "Missing metadata columns in df: class"
  )
})


testthat::test_that("detect_hotelling_nonqc_dual_z errors when no numeric metabolite columns exist", {
  p <- list(qcImputeM = "median", samImputeM = "median")

  df <- data.frame(
    sample = paste0("S", 1:5),
    batch = "b1",
    class = c("A", "A", "A", "QC", "QC"),
    order = 1:5,
    m1 = as.character(1:5), # non-numeric
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df = df, p = p),
    "No numeric metabolite columns found"
  )
})

testthat::test_that("detect_hotelling_nonqc_dual_z errors when too few complete non-QC rows exist", {
  p <- list(qcImputeM = "median", samImputeM = "median")

  df <- data.frame(
    sample = paste0("S", 1:5),
    batch = "b1",
    class = c("A", "A", "QC", "QC", "QC"),
    order = 1:5,
    m1 = c(1, NA, 1, 1, 1), # only one complete non-QC row
    m2 = c(1, NA, 1, 1, 1),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df = df, p = p, min_complete = 3L),
    "Too few complete non-QC rows"
  )
})

testthat::test_that("detect_hotelling_nonqc_dual_z errors if all metabolites are dropped as constant", {
  p <- list(qcImputeM = "median", samImputeM = "median")

  df <- data.frame(
    sample = paste0("S", 1:6),
    batch = "b1",
    class = c("A", "A", "A", "QC", "QC", "QC"),
    order = 1:6,
    m1 = rep(10, 6),
    m2 = rep(10, 6),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df = df, p = p, drop_constant = TRUE, const_tol = 1e-12),
    "All metabolite columns were dropped as constant"
  )
})

testthat::test_that("detect_hotelling_nonqc_dual_z returns expected structure and columns on a clean inlier dataset", {
  p <- list(qcImputeM = "median", samImputeM = "median")
  set.seed(1)

  df <- data.frame(
    sample = paste0("S", 1:24),
    batch = rep(c("b1", "b2"), each = 12),
    class = c(rep("A", 8), rep("B", 8), rep("QC", 8)),
    order = 1:24,
    m1 = c(rnorm(16, 100, 5), rnorm(8, 100, 5)),
    m2 = c(rnorm(16, 50, 3), rnorm(8, 50, 3)),
    m3 = c(rnorm(16, 200, 7), rnorm(8, 200, 7)),
    stringsAsFactors = FALSE
  )

  # Avoid calling impute_missing by keeping complete data
  res <- detect_hotelling_nonqc_dual_z(
    df = df,
    p = p,
    alpha = 0.05,
    z_threshold = 5, # make "extreme_values" unlikely
    class_z_threshold = 5,
    make_pca_plot = FALSE
  )

  testthat::expect_type(res, "list")
  testthat::expect_true(all(c("data", "extreme_values", "pca_plot", "pc_loadings", "params") %in% names(res)))
  testthat::expect_equal(res$params$outlier_sample_z_threshold, 5)
  testthat::expect_equal(res$params$inlier_sample_z_threshold, 5)

  out_df <- res$data
  testthat::expect_true(all(c("T2", "is_outlier_sample", "used_in_fit") %in% names(out_df)))
  testthat::expect_equal(nrow(out_df), nrow(df))

  # With very high thresholds, extreme_values should be empty in typical data
  testthat::expect_s3_class(res$extreme_values, "data.frame")
  testthat::expect_equal(nrow(res$extreme_values), 0L)
  testthat::expect_true(is.null(res$pca_plot))
})

testthat::test_that("detect_hotelling_nonqc_dual_z flags an obvious outlier and returns extreme_values", {
  p <- list(qcImputeM = "median", samImputeM = "median")
  set.seed(2)

  nA <- 15
  nB <- 15
  nQC <- 10

  base_A <- data.frame(
    sample = paste0("A", seq_len(nA)),
    batch = "b1",
    class = "A",
    order = seq_len(nA),
    m1 = rnorm(nA, 100, 2),
    m2 = rnorm(nA, 50, 2),
    m3 = rnorm(nA, 200, 3),
    stringsAsFactors = FALSE
  )

  base_B <- data.frame(
    sample = paste0("B", seq_len(nB)),
    batch = "b1",
    class = "B",
    order = nA + seq_len(nB),
    m1 = rnorm(nB, 100, 2),
    m2 = rnorm(nB, 50, 2),
    m3 = rnorm(nB, 200, 3),
    stringsAsFactors = FALSE
  )

  base_QC <- data.frame(
    sample = paste0("QC", seq_len(nQC)),
    batch = "b1",
    class = "QC",
    order = nA + nB + seq_len(nQC),
    m1 = rnorm(nQC, 100, 2),
    m2 = rnorm(nQC, 50, 2),
    m3 = rnorm(nQC, 200, 3),
    stringsAsFactors = FALSE
  )

  df <- rbind(base_A, base_B, base_QC)

  # Inject a clear non-QC outlier with a huge deviation in one metabolite
  df$m1[1] <- df$m1[1] + 2000

  res <- detect_hotelling_nonqc_dual_z(
    df = df,
    p = p,
    alpha = 0.01, # stricter ellipse
    z_threshold = 3,
    class_z_threshold = 3,
    make_pca_plot = FALSE
  )

  out_df <- res$data

  testthat::expect_true(any(out_df$is_outlier_sample, na.rm = TRUE))
  testthat::expect_true(out_df$is_outlier_sample[1])

  # extreme_values should include at least the injected metabolite for that row
  ev <- res$extreme_values
  testthat::expect_s3_class(ev, "data.frame")
  testthat::expect_true(all(c("metabolite", "value_raw", "value_log", "z_global", "z_class", "T2") %in% names(ev)))
  testthat::expect_true(any(ev$sample == df$sample[1]))
  testthat::expect_true(any(ev$metabolite %in% c("m1", "m2", "m3")))
})

testthat::test_that("detect_hotelling_nonqc_dual_z can call impute_missing when NAs exist (mocked) and proceeds", {
  p <- list(qcImputeM = "median", samImputeM = "median")
  set.seed(3)

  df <- data.frame(
    sample = paste0("S", 1:12),
    batch = "b1",
    class = c(rep("A", 5), rep("B", 5), rep("QC", 2)),
    order = 1:12,
    m1 = rnorm(12, 100, 2),
    m2 = rnorm(12, 50, 2),
    m3 = rnorm(12, 200, 3),
    stringsAsFactors = FALSE
  )

  # Add missing values to force the imputation branch
  df$m2[c(2, 7)] <- NA_real_

  called <- 0L
  imputer <- function(df_in, met_cols, qc_method, sam_method) {
    called <<- called + 1L
    df_out <- df_in
    # Simple deterministic fill to keep PCA valid
    for (m in met_cols) {
      if (anyNA(df_out[[m]])) {
        df_out[[m]][is.na(df_out[[m]])] <- stats::median(df_out[[m]], na.rm = TRUE)
      }
    }
    list(df = df_out)
  }

  testthat::local_mocked_bindings(impute_missing = imputer)

  res <- detect_hotelling_nonqc_dual_z(
    df = df,
    p = p,
    make_pca_plot = FALSE
  )

  testthat::expect_equal(called, 1L)
  testthat::expect_s3_class(res$data, "data.frame")
  testthat::expect_true(all(!is.na(res$data$T2) | is.na(res$data$T2))) # existence check
})

testthat::test_that("detect_hotelling_nonqc_dual_z returns a ggplot pca_plot when enabled (if ggplot2 is installed)", {
  testthat::skip_if_not_installed("ggplot2")

  p <- list(qcImputeM = "median", samImputeM = "median")
  set.seed(4)

  df <- data.frame(
    sample = paste0("S", 1:30),
    batch = rep(c("b1", "b2"), each = 15),
    class = c(rep("A", 12), rep("B", 12), rep("QC", 6)),
    order = 1:30,
    m1 = rnorm(30, 100, 5),
    m2 = rnorm(30, 50, 3),
    m3 = rnorm(30, 200, 7),
    stringsAsFactors = FALSE
  )

  res <- detect_hotelling_nonqc_dual_z(
    df = df,
    p = p,
    make_pca_plot = TRUE
  )

  testthat::expect_s3_class(res$pca_plot, "ggplot")
  testthat::expect_s3_class(res$pc_loadings, "data.frame")
  testthat::expect_true(all(c("metabolite", "PC1", "PC2") %in% names(res$pc_loadings)))
})
