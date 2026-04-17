#' Hotelling T^2 on pooled non-QC samples with global z-score candidate flags
#'
#' This function:
#' 1) extracts metabolite columns,
#' 2) log-transforms and scales them using pooled non-QC samples,
#' 3) fits a PCA model on pooled non-QC samples with complete data,
#' 4) computes Hotelling-like T^2 (via 2D Mahalanobis distance in PC1–PC2
#'    space) for ALL samples (QC and non-QC) with complete data, using the
#'    non-QC PCA model,
#' 5) flags outlier samples based on a chi-square cutoff (df = 2),
#' 6) computes:
#'    - global z-scores (pooled non-QC scaling),
#'    - class-based z-scores (within each non-QC class),
#' 7) builds a candidate outlier table using:
#'    - |global z| >= outlier_sample_z_threshold for samples outside the ellipse,
#'    - |global z| >= inlier_sample_z_threshold for all other non-QC samples,
#'    - QC samples excluded from the candidate table,
#' 8) optionally returns a PCA score plot (PC1 vs PC2) with a 95% Hotelling-like
#'    ellipse drawn on pooled non-QC scores, and
#' 9) returns PC loadings for PC1 and PC2.
#'
#' @param df data.frame
#'   Input data with columns for metadata and metabolite intensities.
#' @param p list
#'   Parameter list used by imputation.
#' @param meta_cols character
#'   Names of metadata columns to exclude from multivariate analysis.
#'   Default: c("sample", "batch", "class", "order").
#' @param class_col character
#'   Name of the column containing class labels. Default: "class".
#' @param qc_label character
#'   Label in `class_col` indicating QC samples. Default: "QC".
#'   QCs are excluded from fitting but included when computing T^2.
#' @param alpha numeric
#'   Significance level for the Hotelling ellipse (default 0.05 => 95% ellipse).
#' @param log_transform logical
#'   Whether to apply log2(x + offset) to metabolite columns. Default TRUE.
#' @param log_offset numeric
#'   Constant added before log-transform to avoid log(0). Default 1.
#' @param outlier_sample_z_threshold numeric
#'   Absolute global z-score threshold for samples outside the Hotelling ellipse.
#'   Default 3.
#' @param inlier_sample_z_threshold numeric
#'   Absolute global z-score threshold for non-QC samples inside the Hotelling
#'   ellipse. Default 5.
#' @param min_complete integer
#'   Minimum number of complete non-QC rows required for PCA / covariance
#'   estimation. Default 3.
#' @param drop_constant logical
#'   If TRUE, drop metabolite columns with near-zero variance (based on
#'   complete non-QC rows) before PCA / covariance estimation. Default TRUE.
#' @param const_tol numeric
#'   Variance threshold below which a metabolite is considered constant.
#'   Default 1e-12.
#' @param ridge_factor numeric
#'   Factor controlling the ridge size added to the covariance matrix in
#'   2D PC space: ridge = ridge_factor * mean(diag(cov_scores)). Default 1e-6.
#' @param make_pca_plot logical
#'   If TRUE (default), compute and return a ggplot object with the PCA scores
#'   and a 95% Hotelling-like ellipse on non-QC samples.
#'
#' @return list
#'   A list with components:
#'   - data: original df with added columns:
#'       * T2: squared Mahalanobis distance in PC1–PC2 space for all samples
#'             with complete data (QC and non-QC; NA otherwise)
#'       * is_outlier_sample: TRUE if sample outside (1 - alpha) ellipse
#'       * used_in_fit: TRUE if row used to fit the PCA / covariance (non-QC)
#'   - extreme_values: long-format data.frame of candidate metabolite values
#'       flagged using global z-score rules
#'   - z_global: wide data.frame of pooled non-QC z-scores
#'   - z_class: wide data.frame of within-class z-scores
#'   - pca_plot: ggplot object with PC1 vs PC2 and a 95% ellipse (or NULL)
#'   - pc_loadings: data.frame with PC1 and PC2 loadings for each retained
#'       metabolite (or NULL if PCA could not be computed)
#'   - params: list of settings used
#'
#' @keywords internal
#' @noRd
detect_hotelling_nonqc_dual_z <- function(
    df,
    p,
    meta_cols                  = c("sample", "batch", "class", "order"),
    class_col                  = "class",
    qc_label                   = "QC",
    alpha                      = 0.05,
    log_transform              = TRUE,
    log_offset                 = 1,
    outlier_sample_z_threshold = 3,
    inlier_sample_z_threshold  = 5,
    min_complete               = 3L,
    drop_constant              = TRUE,
    const_tol                  = 1e-12,
    ridge_factor               = 1e-6,
    make_pca_plot              = TRUE
) {
  missing_meta <- setdiff(meta_cols, names(df))
  if (length(missing_meta) > 0L) {
    stop("Missing metadata columns in df: ",
         paste(missing_meta, collapse = ", "))
  }
  if (!class_col %in% names(df)) {
    stop("class_col '", class_col, "' not found in df.")
  }
  
  candidate_cols <- setdiff(names(df), meta_cols)
  met_cols <- candidate_cols[vapply(df[candidate_cols], is.numeric, logical(1))]
  if (length(met_cols) == 0L) {
    stop("No numeric metabolite columns found.")
  }
  
  if (any(is.na(df[, met_cols, drop = FALSE]))) {
    results <- impute_missing(df, met_cols, p$qcImputeM, p$samImputeM)
    df <- results$df
  }
  
  X_raw <- as.matrix(df[, met_cols, drop = FALSE])
  
  if (log_transform) {
    X_log <- log2(X_raw + log_offset)
  } else {
    X_log <- X_raw
  }
  
  class_vec <- df[[class_col]]
  nonqc_mask <- !is.na(class_vec) & class_vec != qc_label
  
  complete_nonqc <- nonqc_mask & stats::complete.cases(X_log)
  X_fit <- X_log[complete_nonqc, , drop = FALSE]
  
  if (sum(complete_nonqc) < min_complete) {
    stop("Too few complete non-QC rows to estimate PCA / covariance.")
  }
  
  if (drop_constant) {
    v <- apply(X_fit, 2L, stats::var, na.rm = TRUE)
    const_mask <- v <= const_tol | is.na(v)
    
    if (any(const_mask)) {
      dropped_mets <- met_cols[const_mask]
      message(
        "Dropping ", length(dropped_mets),
        " metabolite(s) with near-zero variance among non-QC: ",
        paste(dropped_mets, collapse = ", ")
      )
      
      keep_mask <- !const_mask
      met_cols  <- met_cols[keep_mask]
      X_raw     <- X_raw[, keep_mask, drop = FALSE]
      X_log     <- X_log[, keep_mask, drop = FALSE]
      X_fit     <- X_fit[, keep_mask, drop = FALSE]
    }
  }
  
  if (length(met_cols) == 0L) {
    stop("All metabolite columns were dropped as constant; cannot proceed.")
  }
  
  X_fit_scaled <- scale(X_fit)
  center_scaled <- attr(X_fit_scaled, "scaled:center")
  scale_scaled  <- attr(X_fit_scaled, "scaled:scale")
  
  bad_scale <- !is.finite(scale_scaled) | scale_scaled == 0
  scale_scaled[bad_scale] <- NA_real_
  
  scale_all_rows <- function(X_raw, center_raw, scale_raw) {
    sweep(sweep(X_raw, 2L, center_raw, FUN = "-"), 2L, scale_raw, FUN = "/")
  }
  
  X_scaled_all <- scale_all_rows(X_log, center_scaled, scale_scaled)
  
  n <- nrow(df)
  T2 <- rep(NA_real_, n)
  used_in_fit <- rep(FALSE, n)
  
  complete_nonqc_retained <- nonqc_mask & stats::complete.cases(X_log)
  idx_fit <- which(complete_nonqc_retained)
  
  pca <- NULL
  scores_fit <- NULL
  scores_all <- NULL
  center_scores <- NULL
  cov_scores_reg <- NULL
  pc_loadings <- NULL
  
  if (length(idx_fit) >= 3L && ncol(X_scaled_all) >= 2L) {
    X_nonqc_scaled <- X_scaled_all[idx_fit, , drop = FALSE]
    
    pca <- stats::prcomp(X_nonqc_scaled, center = FALSE, scale. = FALSE)
    scores_fit <- pca$x[, 1:2, drop = FALSE]
    
    center_scores <- colMeans(scores_fit)
    cov_scores <- stats::cov(scores_fit)
    
    d_pc <- diag(cov_scores)
    mean_diag_pc <- mean(d_pc)
    if (!is.finite(mean_diag_pc) || mean_diag_pc <= 0) {
      mean_diag_pc <- 1
    }
    ridge_pc <- ridge_factor * mean_diag_pc
    cov_scores_reg <- cov_scores + diag(ridge_pc, nrow(cov_scores))
    
    pc_loadings <- data.frame(
      metabolite = met_cols,
      PC1 = pca$rotation[, 1],
      PC2 = pca$rotation[, 2],
      stringsAsFactors = FALSE
    )
    
    complete_all_retained <- stats::complete.cases(X_log)
    idx_all <- which(complete_all_retained)
    
    X_all_scaled_complete <- X_scaled_all[idx_all, , drop = FALSE]
    scores_all <- X_all_scaled_complete %*% pca$rotation[, 1:2, drop = FALSE]
    
    T2[idx_all] <- stats::mahalanobis(
      x = scores_all,
      center = center_scores,
      cov = cov_scores_reg
    )
    
    used_in_fit[idx_fit] <- TRUE
  }
  
  cutoff <- stats::qchisq(1 - alpha, df = 2L)
  is_outlier_sample <- !is.na(T2) & (T2 > cutoff)
  
  Z_global <- X_scaled_all
  colnames(Z_global) <- met_cols
  
  Z_class <- matrix(NA_real_, nrow = nrow(df), ncol = length(met_cols))
  colnames(Z_class) <- met_cols
  
  nonqc_classes <- sort(unique(class_vec[nonqc_mask]))
  
  for (cls in nonqc_classes) {
    idx_cls <- which(class_vec == cls)
    X_cls <- X_log[idx_cls, , drop = FALSE]
    
    mu_cls <- apply(X_cls, 2L, mean, na.rm = TRUE)
    sd_cls <- apply(X_cls, 2L, stats::sd, na.rm = TRUE)
    
    zero_sd <- !is.finite(sd_cls) | sd_cls == 0
    sd_cls[zero_sd] <- NA_real_
    
    Z_cls <- sweep(X_cls, 2L, mu_cls, FUN = "-")
    Z_cls <- sweep(Z_cls, 2L, sd_cls, FUN = "/")
    
    Z_class[idx_cls, ] <- Z_cls
  }
  
  # Candidate metabolite values:
  # - outlier non-QC samples use threshold 3
  # - inlier non-QC samples use threshold 5
  # - QC samples excluded
  candidate_sample_mask <- nonqc_mask & !is.na(T2)
  
  candidate_threshold_by_row <- rep(NA_real_, n)
  candidate_threshold_by_row[candidate_sample_mask & is_outlier_sample] <- outlier_sample_z_threshold
  candidate_threshold_by_row[candidate_sample_mask & !is_outlier_sample] <- inlier_sample_z_threshold
  
  candidate_idx <- which(!is.na(candidate_threshold_by_row))
  
  if (length(candidate_idx) > 0L) {
    Zg_candidate <- Z_global[candidate_idx, , drop = FALSE]
    row_thresholds <- matrix(
      candidate_threshold_by_row[candidate_idx],
      nrow = length(candidate_idx),
      ncol = ncol(Zg_candidate)
    )
    
    mask <- !is.na(Zg_candidate) & (abs(Zg_candidate) >= row_thresholds)
    idx <- which(mask, arr.ind = TRUE)
    
    if (nrow(idx) > 0L) {
      row_ids_global <- candidate_idx[idx[, "row"]]
      col_ids <- idx[, "col"]
      metabolite_names <- met_cols[col_ids]
      
      extreme_values <- data.frame(
        df[row_ids_global, meta_cols, drop = FALSE],
        metabolite         = metabolite_names,
        value_raw          = X_raw[cbind(row_ids_global, col_ids)],
        value_log          = X_log[cbind(row_ids_global, col_ids)],
        z_global           = Z_global[cbind(row_ids_global, col_ids)],
        abs_z_global       = abs(Z_global[cbind(row_ids_global, col_ids)]),
        z_class            = Z_class[cbind(row_ids_global, col_ids)],
        abs_z_class        = abs(Z_class[cbind(row_ids_global, col_ids)]),
        T2                 = T2[row_ids_global],
        is_outlier_sample  = is_outlier_sample[row_ids_global],
        z_threshold_used   = candidate_threshold_by_row[row_ids_global],
        stringsAsFactors   = FALSE
      )
      
      extreme_values <- extreme_values[
        order(
          -extreme_values$is_outlier_sample,
          -extreme_values$abs_z_global,
          -extreme_values$T2
        ),
        ,
        drop = FALSE
      ]
      rownames(extreme_values) <- NULL
    } else {
      extreme_values <- data.frame(
        df[0, meta_cols, drop = FALSE],
        metabolite         = character(0),
        value_raw          = numeric(0),
        value_log          = numeric(0),
        z_global           = numeric(0),
        abs_z_global       = numeric(0),
        z_class            = numeric(0),
        abs_z_class        = numeric(0),
        T2                 = numeric(0),
        is_outlier_sample  = logical(0),
        z_threshold_used   = numeric(0),
        stringsAsFactors   = FALSE
      )
    }
  } else {
    extreme_values <- data.frame(
      df[0, meta_cols, drop = FALSE],
      metabolite         = character(0),
      value_raw          = numeric(0),
      value_log          = numeric(0),
      z_global           = numeric(0),
      abs_z_global       = numeric(0),
      z_class            = numeric(0),
      abs_z_class        = numeric(0),
      T2                 = numeric(0),
      is_outlier_sample  = logical(0),
      z_threshold_used   = numeric(0),
      stringsAsFactors   = FALSE
    )
  }
  
  out_df <- df
  out_df$T2 <- T2
  out_df$is_outlier_sample <- is_outlier_sample
  out_df$used_in_fit <- used_in_fit
  
  z_global_df <- data.frame(
    df[, meta_cols, drop = FALSE],
    as.data.frame(Z_global, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  z_class_df <- data.frame(
    df[, meta_cols, drop = FALSE],
    as.data.frame(Z_class, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  pca_plot <- NULL
  
  big_font_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title   = ggplot2::element_text(size = 14, face = "bold"),
      axis.text    = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text  = ggplot2::element_text(size = 10)
    )
  
  .make_hotelling_ellipse <- function(center, cov_mat, cutoff, n = 200L) {
    if (!is.numeric(center) || length(center) != 2L) {
      stop("`center` must be a numeric vector of length 2.")
    }
    if (!is.matrix(cov_mat) || !all(dim(cov_mat) == c(2L, 2L))) {
      stop("`cov_mat` must be a 2x2 numeric matrix.")
    }
    if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0) {
      stop("`cutoff` must be a single positive finite number.")
    }
    if (!is.numeric(n) || length(n) != 1L || n < 4L) {
      stop("`n` must be an integer >= 4.")
    }
    
    eig <- eigen(cov_mat, symmetric = TRUE)
    
    if (any(!is.finite(eig$values)) || any(eig$values <= 0)) {
      stop("Covariance matrix for ellipse must be positive definite.")
    }
    
    theta <- seq(0, 2 * pi, length.out = n)
    unit_circle <- rbind(cos(theta), sin(theta))
    
    transform <- eig$vectors %*% diag(sqrt(eig$values * cutoff), nrow = 2L)
    pts <- sweep(transform %*% unit_circle, 1L, center, FUN = "+")
    
    data.frame(
      PC1 = pts[1L, ],
      PC2 = pts[2L, ],
      stringsAsFactors = FALSE
    )
  }
  
  if (make_pca_plot && !is.null(pca)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' is required for PCA plotting but is not installed.")
    } else {
      PC_mat <- matrix(NA_real_, nrow = n, ncol = 2L)
      colnames(PC_mat) <- c("PC1", "PC2")
      
      complete_all_retained <- stats::complete.cases(X_log)
      idx_all <- which(complete_all_retained)
      
      if (!is.null(scores_all)) {
        PC_mat[idx_all, ] <- scores_all
      }
      
      plot_df <- data.frame(
        sample   = df[[meta_cols[1L]]],
        class    = class_vec,
        is_qc    = class_vec == qc_label,
        is_out   = is_outlier_sample,
        PC1      = PC_mat[, "PC1"],
        PC2      = PC_mat[, "PC2"],
        stringsAsFactors = FALSE
      )
      
      plot_df$group <- "Non-QC inlier"
      plot_df$group[plot_df$is_qc] <- "QC"
      plot_df$group[plot_df$is_out] <- "Outside ellipse"
      plot_df$group <- factor(
        plot_df$group,
        levels = c("QC", "Non-QC inlier", "Outside ellipse")
      )
      
      ellipse_df <- .make_hotelling_ellipse(
        center = center_scores,
        cov_mat = cov_scores_reg,
        cutoff = cutoff,
        n = 300L
      )
      
      pca_plot <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = PC1, y = PC2)
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = group),
          alpha = 0.8,
          size = 2
        ) +
        ggplot2::geom_path(
          data = ellipse_df,
          mapping = ggplot2::aes(x = PC1, y = PC2),
          inherit.aes = FALSE,
          linetype = "dashed"
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "QC" = "#999999",
            "Non-QC inlier" = "#1f78b4",
            "Outside ellipse" = "#e31a1c"
          )
        ) +
        big_font_theme +
        ggplot2::labs(
          title = "PCA (PC1–PC2) with Mahalanobis 95% ellipse (fit on non-QC)",
          color = "Group"
        )
    }
  }
  
  list(
    data           = out_df,
    extreme_values = extreme_values,
    z_global       = z_global_df,
    z_class        = z_class_df,
    pca_plot       = pca_plot,
    pc_loadings    = pc_loadings,
    params         = list(
      alpha                      = alpha,
      cutoff                     = cutoff,
      outlier_sample_z_threshold = outlier_sample_z_threshold,
      inlier_sample_z_threshold  = inlier_sample_z_threshold,
      log_transform              = log_transform,
      log_offset                 = log_offset,
      p                          = ncol(X_scaled_all),
      class_col                  = class_col,
      qc_label                   = qc_label,
      min_complete               = min_complete,
      drop_constant              = drop_constant,
      const_tol                  = const_tol,
      ridge_factor               = ridge_factor,
      retained_metabolites       = met_cols
    )
  )
}