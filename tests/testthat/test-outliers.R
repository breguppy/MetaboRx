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
    m2 = 50 + rep((i %% 4L) * 2, 3L) +
      rep(c(0, 1, -4), each = n),
    m3 = 200 + rep((i %% 5L) * 3, 3L) +
      rep(c(0, -2, 6), each = n),
    stringsAsFactors = FALSE
  )
}

# Each metabolite has exactly one qualifying value in each scoring pool.
# Identical columns also provide exact ties for testing selection order.
hotelling_selection_data <- function(n_metabolites = 20L) {
  n <- 40L
  
  df <- data.frame(
    sample = paste0("S", seq_len(2L * n)),
    batch = "b1",
    class = c(rep(c("A", "B"), each = n / 2L), rep("QC", n)),
    order = seq_len(2L * n),
    stringsAsFactors = FALSE
  )
  
  study_values <- as.numeric(seq_len(n))
  study_values[1L] <- 10000
  
  qc_values <- as.numeric(seq_len(n))
  qc_values[1L] <- -10000
  
  for (j in seq_len(n_metabolites)) {
    df[[paste0("m", j)]] <- c(study_values, qc_values)
  }
  
  df
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
      df,
      hotelling_test_params(),
      min_complete = 3L
    ),
    "Too few complete non-QC rows"
  )
})

testthat::test_that("the result includes pooled scores and typed empty candidates", {
  df <- hotelling_test_data()
  res <- detect_hotelling_nonqc_dual_z(
    df,
    hotelling_test_params(),
    make_pca_plot = FALSE
  )
  
  testthat::expect_named(
    res,
    c(
      "data", "extreme_values", "z_global", "z_class",
      "pca_plot", "pc_loadings", "params"
    ),
    ignore.order = TRUE
  )
  testthat::expect_equal(res$params$value_z_threshold, 2)
  testthat::expect_equal(res$params$flag_fraction, 0.10)
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
      "z_threshold_used", "flag_fraction"
    ),
    ignore.order = TRUE
  )
  testthat::expect_type(ev$z_global, "double")
  testthat::expect_type(ev$T2, "double")
  testthat::expect_type(ev$is_outlier_sample, "logical")
})

testthat::test_that("study samples are pooled and QCs are scored separately", {
  df <- hotelling_test_data()
  res <- detect_hotelling_nonqc_dual_z(
    df,
    hotelling_test_params(),
    make_pca_plot = FALSE
  )
  
  groups <- list(
    study = which(df$class != "QC"),
    qc = which(df$class == "QC")
  )
  
  for (rows in groups) {
    for (metabolite in c("m1", "m2", "m3")) {
      expected <- as.numeric(scale(log2(df[[metabolite]][rows] + 1)))
      testthat::expect_equal(
        res$z_global[[metabolite]][rows],
        expected
      )
    }
  }
})

testthat::test_that("PCA uses study rows and distances use study and QC scores", {
  df <- hotelling_test_data()
  study <- df$class != "QC"
  qc_row <- which(!study)[1L]
  df$m1[qc_row] <- 500
  
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

testthat::test_that("the largest 10 percent above 2 are selected separately per pool", {
  df <- hotelling_selection_data()
  
  res <- detect_hotelling_nonqc_dual_z(
    df,
    hotelling_test_params(),
    log_transform = FALSE,
    alpha = 1e-12,
    make_pca_plot = FALSE
  )
  
  metabolites <- paste0("m", seq_len(20L))
  study <- df$class != "QC"
  
  # There are 20 qualifying values per pool, so two must be selected.
  study_z <- as.matrix(res$z_global[study, metabolites])
  qc_z <- as.matrix(res$z_global[!study, metabolites])
  
  testthat::expect_equal(sum(abs(study_z) > 2), 20L)
  testthat::expect_equal(sum(abs(qc_z) > 2), 20L)
  
  ev <- res$extreme_values
  study_ev <- ev[ev$class != "QC", , drop = FALSE]
  qc_ev <- ev[ev$class == "QC", , drop = FALSE]
  
  testthat::expect_equal(nrow(study_ev), 2L)
  testthat::expect_equal(nrow(qc_ev), 2L)
  
  # Exact ties resolve by input row, then input metabolite column.
  testthat::expect_setequal(study_ev$metabolite, c("m1", "m2"))
  testthat::expect_setequal(qc_ev$metabolite, c("m1", "m2"))
  testthat::expect_true(all(study_ev$sample == "S1"))
  testthat::expect_true(all(qc_ev$sample == "S41"))
  
  testthat::expect_true(all(study_ev$z_global > 2))
  testthat::expect_true(all(qc_ev$z_global < -2))
  testthat::expect_equal(ev$abs_z_global, abs(ev$z_global))
  testthat::expect_equal(ev$z_class, ev$z_global)
  testthat::expect_equal(ev$abs_z_class, ev$abs_z_global)
  testthat::expect_equal(ev$z_threshold_used, rep(2, nrow(ev)))
  testthat::expect_equal(ev$flag_fraction, rep(0.10, nrow(ev)))
  
  # The wide ellipse keeps these candidates inside the sample cutoff.
  testthat::expect_false(any(ev$is_outlier_sample))
})

testthat::test_that("selection counts round down after thresholding", {
  for (n_metabolites in c(9L, 19L)) {
    df <- hotelling_selection_data(n_metabolites)
    
    res <- detect_hotelling_nonqc_dual_z(
      df,
      hotelling_test_params(),
      log_transform = FALSE,
      make_pca_plot = FALSE
    )
    
    ev <- res$extreme_values
    expected_per_pool <- floor(n_metabolites * 0.10)
    
    testthat::expect_equal(
      sum(ev$class != "QC"),
      expected_per_pool
    )
    testthat::expect_equal(
      sum(ev$class == "QC"),
      expected_per_pool
    )
  }
})

testthat::test_that("a metabolite constant within QC has missing QC Z-scores", {
  df <- hotelling_test_data()
  qc <- df$class == "QC"
  df$m1[qc] <- 110
  
  res <- detect_hotelling_nonqc_dual_z(
    df,
    hotelling_test_params(),
    make_pca_plot = FALSE
  )
  
  testthat::expect_true(all(is.na(res$z_global$m1[qc])))
  testthat::expect_true(all(is.finite(res$z_global$m1[!qc])))
  testthat::expect_false(
    any(
      res$extreme_values$class == "QC" &
        res$extreme_values$metabolite == "m1"
    )
  )
})

testthat::test_that("metabolites excluded from PCA still receive pooled Z-scores", {
  df <- hotelling_test_data()
  study <- df$class != "QC"
  df[study, c("m1", "m2", "m3")] <- 10
  
  res <- suppressMessages(
    detect_hotelling_nonqc_dual_z(
      df,
      hotelling_test_params(),
      make_pca_plot = FALSE
    )
  )
  
  testthat::expect_null(res$pc_loadings)
  testthat::expect_null(res$pca_plot)
  testthat::expect_true(all(is.na(res$data$T2)))
  testthat::expect_false(any(res$data$used_in_fit))
  testthat::expect_true(all(is.na(res$z_global$m1[study])))
  testthat::expect_true(all(is.finite(res$z_global$m1[!study])))
  testthat::expect_setequal(
    res$params$screened_metabolites,
    c("m1", "m2", "m3")
  )
  testthat::expect_length(res$params$retained_metabolites, 0L)
})

testthat::test_that("imputation supports PCA but observed values determine Z-scores", {
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
  
  res <- detect_hotelling_nonqc_dual_z(
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
  testthat::expect_true(all(is.na(res$z_global$m2[missing_rows])))
  
  for (rows in list(which(df$class != "QC"), which(df$class == "QC"))) {
    observed_rows <- rows[!is.na(df$m2[rows])]
    expected <- as.numeric(
      scale(log2(df$m2[observed_rows] + 1))
    )
    
    testthat::expect_equal(
      res$z_global$m2[observed_rows],
      expected
    )
  }
  
  ev <- res$extreme_values
  testthat::expect_false(
    any(ev$sample %in% df$sample[missing_rows] & ev$metabolite == "m2")
  )
})

testthat::test_that("the PCA plot and loadings are returned when enabled", {
  testthat::skip_if_not_installed("ggplot2")
  df <- hotelling_test_data()
  
  res <- detect_hotelling_nonqc_dual_z(
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