# tests/testthat/test-detect_outliers.R

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
    m2 = 50 + rep((i %% 4L) * 2, 3L) +
      rep(c(0, 1, -4), each = n),
    m3 = 200 + rep((i %% 5L) * 3, 3L) +
      rep(c(0, -2, 6), each = n),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("required metadata columns must be present", {
  df <- hotelling_test_data()
  p <- hotelling_test_params()
  
  df$order <- NULL
  testthat::expect_error(
    detect_outliers(df, p),
    "Missing metadata columns in df: order"
  )
  
  df$order <- seq_len(nrow(df))
  df$class <- NULL
  testthat::expect_error(
    detect_outliers(df, p),
    "Missing metadata columns in df: class"
  )
})

testthat::test_that("at least one numeric metabolite column is required", {
  df <- hotelling_test_data()
  df$m1 <- as.character(df$m1)
  df$m2 <- NULL
  df$m3 <- NULL
  
  testthat::expect_error(
    detect_outliers(df, hotelling_test_params()),
    "No numeric metabolite columns found"
  )
})

testthat::test_that("QC rows cannot satisfy the minimum study row count", {
  df <- hotelling_test_data()
  df <- df[c(1L, 2L, which(df$class == "QC")), , drop = FALSE]
  
  testthat::expect_error(
    detect_outliers(
      df,
      hotelling_test_params(),
      min_complete = 3L
    ),
    "Too few complete non-QC rows"
  )
})

testthat::test_that("the result includes method-specific scores and typed empty candidates", {
  df <- hotelling_test_data()
  res <- detect_outliers(
    df,
    hotelling_test_params(),
    make_pca_plot = FALSE
  )
  
  testthat::expect_named(
    res,
    c(
      "data", "extreme_values", "rout_scores", "value_scores",
      "z_global", "z_class",
      "pca_plot", "pc_loadings", "params"
    ),
    ignore.order = TRUE
  )
  testthat::expect_equal(res$params$grubbs_alpha, 0.005)
  testthat::expect_equal(res$params$rout_q, 0.001)
  testthat::expect_equal(res$params$min_rout_observations, 8L)
  testthat::expect_equal(nrow(res$data), nrow(df))
  testthat::expect_null(res$pca_plot)
  testthat::expect_identical(res$z_class, res$z_global)
  
  ev <- res$extreme_values
  testthat::expect_s3_class(ev, "data.frame")
  testthat::expect_equal(nrow(ev), 0L)
  testthat::expect_named(
    ev,
    c(
      "sample", "batch", "class", "order", "metabolite",
      "value_raw", "value_log", "z_global", "abs_z_global",
      "z_class", "abs_z_class", "T2", "is_outlier_sample",
      "method", "test_statistic", "adjusted_p", "p_cutoff",
      "rout_score", "rout_adjusted_p", "z_threshold_used", "flag_fraction"
    ),
    ignore.order = TRUE
  )
  testthat::expect_type(ev$method, "character")
  testthat::expect_type(ev$test_statistic, "double")
  testthat::expect_type(ev$adjusted_p, "double")
  testthat::expect_type(ev$T2, "double")
  testthat::expect_type(ev$is_outlier_sample, "logical")
})

testthat::test_that("PCA uses study rows and distances use study and QC scores", {
  df <- hotelling_test_data()
  study <- df$class != "QC"
  qc_row <- which(!study)[1L]
  df$m1[qc_row] <- 500
  
  alpha <- 0.05
  ridge_factor <- 1e-6
  
  res <- detect_outliers(
    df,
    hotelling_test_params(),
    alpha = alpha,
    ridge_factor = ridge_factor,
    make_pca_plot = FALSE
  )
  
  x <- log2(as.matrix(df[c("m1", "m2", "m3")]) + 1)
  pca <- stats::prcomp(
    x[study, , drop = FALSE],
    center = TRUE,
    scale. = TRUE
  )
  scores <- stats::predict(pca, newdata = x)[, 1:2, drop = FALSE]
  score_cov <- stats::cov(scores)
  ridge <- ridge_factor * mean(diag(score_cov))
  
  expected_t2 <- stats::mahalanobis(
    scores,
    colMeans(scores),
    score_cov + diag(ridge, 2L)
  )
  cutoff <- stats::qchisq(1 - alpha, df = 2L)
  
  testthat::expect_identical(res$data$used_in_fit, study)
  testthat::expect_equal(res$data$T2, unname(expected_t2))
  testthat::expect_equal(res$params$cutoff, cutoff)
  testthat::expect_identical(
    res$data$is_outlier_sample,
    unname(expected_t2 > cutoff)
  )
  testthat::expect_true(res$data$is_outlier_sample[qc_row])
})

testthat::test_that("seven observed values use Grubbs-type and eight use ROUT", {
  df <- hotelling_test_data()
  rows <- c(which(df$class == "A")[1:8], which(df$class == "B")[1:8],
            which(df$class == "QC"))
  df <- df[rows, , drop = FALSE]
  study_a <- which(df$class == "A")
  study_b <- which(df$class == "B")
  qc <- which(df$class == "QC")
  df$m1 <- 2^(10 + c(
    seq(-0.03, 0.03, length.out = 6L), 5, 0,
    seq(-0.03, 0.03, length.out = 7L), -5,
    rep(0, length(qc))
  )) - 1
  observed <- df
  observed$m1[study_a[8L]] <- NA_real_
  
  res <- detect_outliers(
    df,
    hotelling_test_params(),
    log_transform = FALSE,
    make_pca_plot = FALSE,
    observed_df = observed
  )
  ev <- res$extreme_values
  ev <- ev[ev$metabolite == "m1", , drop = FALSE]
  testthat::expect_setequal(ev$sample, df$sample[c(study_a[7L], study_b[8L])])
  
  a <- ev[ev$method == "Grubbs-type", , drop = FALSE]
  b <- ev[ev$method == "ROUT", , drop = FALSE]
  testthat::expect_identical(a$sample, df$sample[study_a[7L]])
  testthat::expect_identical(b$sample, df$sample[study_b[8L]])
  testthat::expect_equal(a$p_cutoff, 0.005)
  testthat::expect_equal(b$p_cutoff, 0.001)
  testthat::expect_true(all(ev$adjusted_p < ev$p_cutoff))
  testthat::expect_equal(ev$value_log, log2(ev$value_raw + 1))
  
  # Check the small-class result against a one-way model using all classes.
  reference <- observed[is.finite(observed$m1), , drop = FALSE]
  fit <- stats::lm(log2(m1 + 1) ~ class, data = reference)
  scores <- stats::rstudent(fit)[reference$class == "A"]
  p_values <- 2 * stats::pt(-abs(scores), df = stats::df.residual(fit) - 1L)
  expected <- stats::p.adjust(p_values, method = "bonferroni")
  testthat::expect_equal(res$value_scores$m1[study_a[1:7]], unname(scores))
  testthat::expect_equal(a$adjusted_p, unname(expected[7L]))
  testthat::expect_lt(b$test_statistic, 0)
  testthat::expect_true(is.na(res$value_scores$m1[study_a[8L]]))
  testthat::expect_true(all(is.na(res$rout_scores$m1[c(study_a, qc)])))
  testthat::expect_equal(res$data$m1, df$m1)
})

testthat::test_that("metabolites excluded from PCA still receive value scores", {
  df <- hotelling_test_data()
  study <- df$class != "QC"
  df$m1[study] <- 100
  df$m1[!study] <- 2^(10 + c(seq(-0.05, 0.05, length.out = 11L), 5)) - 1
  
  res <- suppressMessages(detect_outliers(
    df,
    hotelling_test_params(),
    make_pca_plot = FALSE
  ))
  testthat::expect_false("m1" %in% res$pc_loadings$metabolite)
  testthat::expect_true("m1" %in% res$params$screened_metabolites)
  testthat::expect_true(all(is.na(res$value_scores$m1[study])))
  testthat::expect_true(all(is.finite(res$rout_scores$m1[!study])))
  testthat::expect_true(any(
    res$extreme_values$metabolite == "m1" &
      res$extreme_values$class == "QC"
  ))
})

testthat::test_that("imputation supports PCA but observed values determine candidate scores", {
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
  
  res <- detect_outliers(
    df,
    p,
    make_pca_plot = FALSE
  )
  
  testthat::expect_equal(calls, 1L)
  testthat::expect_equal(
    res$data$m2[missing_rows],
    rep(fill, length(missing_rows))
  )
  testthat::expect_true(all(is.finite(res$data$T2)))
  testthat::expect_true(all(is.na(res$value_scores$m2[missing_rows])))
  
  cached <- detect_outliers(
    res$data[, names(df), drop = FALSE],
    p,
    make_pca_plot = FALSE,
    observed_df = df[nrow(df):1L, , drop = FALSE]
  )
  testthat::expect_equal(calls, 1L)
  testthat::expect_equal(cached$value_scores, res$value_scores)
  testthat::expect_equal(cached$extreme_values, res$extreme_values)
  
  ev <- res$extreme_values
  testthat::expect_false(
    any(ev$sample %in% df$sample[missing_rows] & ev$metabolite == "m2")
  )
})

testthat::test_that("the PCA plot and loadings are returned when enabled", {
  testthat::skip_if_not_installed("ggplot2")
  df <- hotelling_test_data()
  
  res <- detect_outliers(
    df,
    hotelling_test_params(),
    make_pca_plot = TRUE
  )
  
  testthat::expect_s3_class(res$pca_plot, "ggplot")
  testthat::expect_equal(nrow(res$pca_plot$data), nrow(df))
  testthat::expect_s3_class(res$pc_loadings, "data.frame")
  testthat::expect_named(res$pc_loadings, c("metabolite", "PC1", "PC2"))
  testthat::expect_setequal(
    res$pc_loadings$metabolite,
    c("m1", "m2", "m3")
  )
})