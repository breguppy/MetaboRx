# tests/testthat/test-detect_hotelling_nonqc_dual_z.R

hotelling_test_params <- function() {
  list(qcImputeM = "median", samImputeM = "median")
}

hotelling_test_data <- function() {
  n <- 12L
  i <- seq_len(n)
  
  data.frame(
    sample = paste0("S", seq_len(3L * n)),
    batch = "b1",
    class = rep(c("A", "B", "QC"), each = n),
    order = seq_len(3L * n),
    m1 = 100 + rep(i, 3L) + rep(c(0, 2, 8), each = n),
    m2 = 50 + rep((i %% 4L) * 2, 3L) + rep(c(0, 1, -4), each = n),
    m3 = 200 + rep((i %% 5L) * 3, 3L) + rep(c(0, -2, 6), each = n),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("required metadata columns must be present", {
  df <- hotelling_test_data()
  p <- hotelling_test_params()
  
  df$order <- NULL
  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df, p),
    "Missing metadata columns in df: order"
  )
  
  df$order <- seq_len(nrow(df))
  df$class <- NULL
  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df, p),
    "Missing metadata columns in df: class"
  )
})

testthat::test_that("at least one numeric metabolite column is required", {
  df <- hotelling_test_data()
  df$m1 <- as.character(df$m1)
  df$m2 <- NULL
  df$m3 <- NULL
  
  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(df, hotelling_test_params()),
    "No numeric metabolite columns found"
  )
})

testthat::test_that("QC rows cannot satisfy the minimum study row count", {
  df <- hotelling_test_data()
  df <- df[c(1L, 2L, which(df$class == "QC")), , drop = FALSE]
  
  testthat::expect_error(
    detect_hotelling_nonqc_dual_z(
      df, hotelling_test_params(), min_complete = 3L
    ),
    "Too few complete non-QC rows"
  )
})

testthat::test_that("the detector errors if all metabolites are constant in study samples", {
  df <- hotelling_test_data()
  study <- df$class != "QC"
  df[study, c("m1", "m2", "m3")] <- 10
  
  testthat::expect_error(
    suppressMessages(detect_hotelling_nonqc_dual_z(df, hotelling_test_params())),
    "All metabolite columns were dropped as constant"
  )
})

testthat::test_that("the result has class scores and a typed empty candidate table", {
  df <- hotelling_test_data()
  res <- detect_hotelling_nonqc_dual_z(
    df, hotelling_test_params(), make_pca_plot = FALSE
  )
  
  testthat::expect_named(
    res,
    c("data", "extreme_values", "z_class", "pca_plot", "pc_loadings", "params"),
    ignore.order = TRUE
  )
  testthat::expect_equal(res$params$class_z_threshold, 3)
  testthat::expect_equal(nrow(res$data), nrow(df))
  testthat::expect_null(res$pca_plot)
  
  testthat::expect_s3_class(res$extreme_values, "data.frame")
  testthat::expect_equal(nrow(res$extreme_values), 0L)
  testthat::expect_named(
    res$extreme_values,
    c(
      "sample", "batch", "class", "order", "metabolite", "value_raw",
      "value_log", "z_class", "abs_z_class", "T2", "is_outlier_sample",
      "z_threshold_used"
    ),
    ignore.order = TRUE
  )
  testthat::expect_type(res$extreme_values$z_class, "double")
  testthat::expect_type(res$extreme_values$T2, "double")
  testthat::expect_type(res$extreme_values$is_outlier_sample, "logical")
  
  for (class_label in unique(df$class)) {
    rows <- df$class == class_label
    expected <- as.numeric(scale(log2(df$m1[rows] + 1)))
    testthat::expect_equal(res$z_class$m1[rows], expected)
  }
})

testthat::test_that("PCA uses study rows and Mahalanobis distances use study and QC scores", {
  df <- hotelling_test_data()
  study <- df$class != "QC"
  qc <- !study
  df$m1[which(qc)[1L]] <- 500
  
  alpha <- 0.05
  ridge_factor <- 1e-6
  res <- detect_hotelling_nonqc_dual_z(
    df,
    hotelling_test_params(),
    alpha = alpha,
    ridge_factor = ridge_factor,
    make_pca_plot = FALSE
  )
  
  x <- log2(as.matrix(df[c("m1", "m2", "m3")]) + 1)
  pca <- stats::prcomp(x[study, , drop = FALSE], center = TRUE, scale. = TRUE)
  scores <- stats::predict(pca, newdata = x)[, 1:2, drop = FALSE]
  score_cov <- stats::cov(scores)
  ridge <- ridge_factor * mean(diag(score_cov))
  expected_t2 <- stats::mahalanobis(
    scores, colMeans(scores), score_cov + diag(ridge, 2L)
  )
  cutoff <- stats::qchisq(1 - alpha, df = 2L)
  
  testthat::expect_identical(res$data$used_in_fit, study)
  testthat::expect_equal(res$data$T2, unname(expected_t2))
  testthat::expect_equal(res$params$cutoff, cutoff)
  testthat::expect_identical(
    res$data$is_outlier_sample, unname(expected_t2 > cutoff)
  )
  testthat::expect_true(res$data$is_outlier_sample[which(qc)[1L]])
})

testthat::test_that("class Z-scores flag high and low values independently of sample flags", {
  df <- hotelling_test_data()
  study_row <- which(df$class == "A")[1L]
  qc_row <- which(df$class == "QC")[1L]
  df$m1[study_row] <- 500
  df$m2[qc_row] <- 5
  
  res <- detect_hotelling_nonqc_dual_z(
    df,
    hotelling_test_params(),
    # Use a wide ellipse so both candidates remain inside it.
    alpha = 1e-12,
    make_pca_plot = FALSE
  )
  
  ev <- res$extreme_values
  study_hit <- which(ev$sample == df$sample[study_row] & ev$metabolite == "m1")
  qc_hit <- which(ev$sample == df$sample[qc_row] & ev$metabolite == "m2")
  
  testthat::expect_length(study_hit, 1L)
  testthat::expect_length(qc_hit, 1L)
  testthat::expect_gt(ev$z_class[study_hit], 3)
  testthat::expect_lt(ev$z_class[qc_hit], -3)
  testthat::expect_false(any(ev$is_outlier_sample[c(study_hit, qc_hit)]))
  
  z <- as.matrix(res$z_class[c("m1", "m2", "m3")])
  expected <- which(is.finite(z) & abs(z) > 3, arr.ind = TRUE)
  testthat::expect_setequal(
    paste(ev$sample, ev$metabolite),
    paste(df$sample[expected[, "row"]], colnames(z)[expected[, "col"]])
  )
  testthat::expect_equal(ev$abs_z_class, abs(ev$z_class))
  testthat::expect_equal(ev$z_threshold_used, rep(3, nrow(ev)))
})

testthat::test_that("a metabolite constant within QC has missing QC Z-scores", {
  df <- hotelling_test_data()
  qc <- df$class == "QC"
  df$m1[qc] <- 110
  
  res <- detect_hotelling_nonqc_dual_z(
    df, hotelling_test_params(), make_pca_plot = FALSE
  )
  
  testthat::expect_true(all(is.na(res$z_class$m1[qc])))
  testthat::expect_true(all(is.finite(res$z_class$m1[!qc])))
  testthat::expect_false(
    any(res$extreme_values$class == "QC" & res$extreme_values$metabolite == "m1")
  )
})

testthat::test_that("missing values are imputed before scoring", {
  df <- hotelling_test_data()
  p <- list(qcImputeM = "median", samImputeM = "mean")
  missing_rows <- c(2L, which(df$class == "QC")[1L])
  df$m2[missing_rows] <- NA_real_
  fill <- stats::median(df$m2, na.rm = TRUE)
  calls <- 0L
  
  imputer <- function(df_in, met_cols, qc_method, sam_method) {
    calls <<- calls + 1L
    testthat::expect_setequal(met_cols, c("m1", "m2", "m3"))
    testthat::expect_identical(qc_method, p$qcImputeM)
    testthat::expect_identical(sam_method, p$samImputeM)
    df_in$m2[is.na(df_in$m2)] <- fill
    list(df = df_in)
  }
  
  testthat::local_mocked_bindings(impute_missing = imputer)
  
  res <- detect_hotelling_nonqc_dual_z(df, p, make_pca_plot = FALSE)
  
  testthat::expect_equal(calls, 1L)
  testthat::expect_equal(res$data$m2[missing_rows], rep(fill, length(missing_rows)))
  testthat::expect_true(all(is.finite(res$data$T2)))
  testthat::expect_true(all(is.finite(res$z_class$m2[missing_rows])))
})

testthat::test_that("the PCA plot and loadings are returned when plotting is enabled", {
  testthat::skip_if_not_installed("ggplot2")
  df <- hotelling_test_data()
  
  res <- detect_hotelling_nonqc_dual_z(
    df, hotelling_test_params(), make_pca_plot = TRUE
  )
  
  testthat::expect_s3_class(res$pca_plot, "ggplot")
  testthat::expect_equal(nrow(res$pca_plot$data), nrow(df))
  testthat::expect_s3_class(res$pc_loadings, "data.frame")
  testthat::expect_named(res$pc_loadings, c("metabolite", "PC1", "PC2"))
  testthat::expect_setequal(res$pc_loadings$metabolite, c("m1", "m2", "m3"))
})