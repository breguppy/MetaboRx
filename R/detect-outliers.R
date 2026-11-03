#' Hotelling T^2 on pooled non-QC samples with global z-score candidate flags
#'
#' This function:
#' 1) extracts metabolite columns,
#' 2) log-transforms and scales them using pooled non-QC samples,
#' 3) fits a PCA model on pooled non-QC samples with complete data,
#' 4) computes Hotelling-like T^2 (via 2D Mahalanobis distance in PC1-PC2
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
#' @param z_threshold,class_z_threshold Deprecated aliases accepted for
#'   backward compatibility with older tests/callers.
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
#'       * T2: squared Mahalanobis distance in PC1-PC2 space for all samples
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
#' Detect sample and metabolite outliers
#'
#' Fits PCA on scaled, complete non-QC rows and estimates the PC1-PC2
#' Mahalanobis reference from all complete projected rows, including QCs.
#' Flags individual retained metabolite values when their absolute
#' within-class Z-score exceeds 3. QCs are treated as their own class.
#'
#' The four legacy Z-score threshold arguments are accepted for compatibility.
#' `used_in_fit` identifies rows used in the PCA fit.
#'
#' @keywords internal
#' @noRd
detect_hotelling_nonqc_dual_z <- function(
    df,
    p,
    meta_cols = c("sample", "batch", "class", "order"),
    class_col = "class",
    qc_label = "QC",
    alpha = 0.05,
    log_transform = TRUE,
    log_offset = 1,
    outlier_sample_z_threshold = 3,
    inlier_sample_z_threshold = 5,
    z_threshold = NULL,
    class_z_threshold = NULL,
    min_complete = 3L,
    drop_constant = TRUE,
    const_tol = 1e-12,
    ridge_factor = 1e-6,
    make_pca_plot = TRUE
) {
  missing_meta <- setdiff(meta_cols, names(df))
  if (length(missing_meta) > 0L) {
    stop("Missing metadata columns in df: ", paste(missing_meta, collapse = ", "))
  }
  if (!class_col %in% names(df)) {
    stop("class_col '", class_col, "' not found in df.")
  }
  
  candidate_cols <- setdiff(names(df), meta_cols)
  met_cols <- candidate_cols[vapply(df[candidate_cols], is.numeric, logical(1))]
  if (length(met_cols) == 0L) {
    stop("No numeric metabolite columns found.")
  }
  
  if (anyNA(df[, met_cols, drop = FALSE])) {
    results <- impute_missing(df, met_cols, p$qcImputeM, p$samImputeM)
    df <- results$df
  }
  
  X_raw <- as.matrix(df[, met_cols, drop = FALSE])
  X_log <- if (log_transform) log2(X_raw + log_offset) else X_raw
  
  class_vec <- df[[class_col]]
  nonqc_mask <- !is.na(class_vec) & class_vec != qc_label
  complete_nonqc <- nonqc_mask & stats::complete.cases(X_log)
  X_fit <- X_log[complete_nonqc, , drop = FALSE]
  
  if (sum(complete_nonqc) < min_complete) {
    stop("Too few complete non-QC rows to estimate PCA.")
  }
  
  if (drop_constant) {
    v <- apply(X_fit, 2L, stats::var, na.rm = TRUE)
    const_mask <- v <= const_tol | is.na(v)
    
    if (any(const_mask)) {
      message(
        "Dropping ", sum(const_mask),
        " metabolite(s) with near-zero variance among non-QC: ",
        paste(met_cols[const_mask], collapse = ", ")
      )
      
      keep_mask <- !const_mask
      met_cols <- met_cols[keep_mask]
      X_raw <- X_raw[, keep_mask, drop = FALSE]
      X_log <- X_log[, keep_mask, drop = FALSE]
      X_fit <- X_fit[, keep_mask, drop = FALSE]
    }
  }
  
  if (length(met_cols) == 0L) {
    stop("All metabolite columns were dropped as constant; cannot proceed.")
  }
  
  X_fit_scaled <- scale(X_fit)
  center_scaled <- attr(X_fit_scaled, "scaled:center")
  scale_scaled <- attr(X_fit_scaled, "scaled:scale")
  scale_scaled[!is.finite(scale_scaled) | scale_scaled == 0] <- NA_real_
  
  X_scaled_all <- sweep(X_log, 2L, center_scaled, FUN = "-")
  X_scaled_all <- sweep(X_scaled_all, 2L, scale_scaled, FUN = "/")
  
  n <- nrow(df)
  T2 <- rep(NA_real_, n)
  used_in_fit <- rep(FALSE, n)
  cutoff <- stats::qchisq(1 - alpha, df = 2L)
  
  complete_scaled <- rowSums(!is.finite(X_scaled_all)) == 0L
  idx_all <- which(complete_scaled)
  idx_fit <- which(nonqc_mask & complete_scaled)
  
  pca <- NULL
  pc_loadings <- NULL
  
  if (length(idx_fit) >= 3L && ncol(X_scaled_all) >= 2L) {
    pca <- stats::prcomp(
      X_scaled_all[idx_fit, , drop = FALSE],
      center = FALSE,
      scale. = FALSE
    )
    
    rotation <- pca$rotation[, 1:2, drop = FALSE]
    scores_all <- X_scaled_all[idx_all, , drop = FALSE] %*% rotation
    
    # Use all complete projected rows, including QCs, for the ellipse.
    center_scores <- colMeans(scores_all)
    cov_scores <- stats::cov(scores_all)
    mean_diag_pc <- mean(diag(cov_scores))
    if (!is.finite(mean_diag_pc) || mean_diag_pc <= 0) {
      mean_diag_pc <- 1
    }
    
    cov_scores_reg <- cov_scores + diag(ridge_factor * mean_diag_pc, 2L)
    T2[idx_all] <- stats::mahalanobis(scores_all, center_scores, cov_scores_reg)
    used_in_fit[idx_fit] <- TRUE
    
    pc_loadings <- data.frame(
      metabolite = met_cols,
      PC1 = rotation[, 1L],
      PC2 = rotation[, 2L],
      stringsAsFactors = FALSE
    )
  }
  
  is_outlier_sample <- !is.na(T2) & T2 > cutoff
  
  Z_class <- matrix(
    NA_real_,
    nrow = n,
    ncol = length(met_cols),
    dimnames = list(NULL, met_cols)
  )
  
  classes <- sort(unique(class_vec[!is.na(class_vec)]))
  for (cls in classes) {
    idx_cls <- which(class_vec == cls)
    X_cls <- X_log[idx_cls, , drop = FALSE]
    
    mu_cls <- colMeans(X_cls, na.rm = TRUE)
    sd_cls <- apply(X_cls, 2L, stats::sd, na.rm = TRUE)
    sd_cls[!is.finite(sd_cls) | sd_cls == 0] <- NA_real_
    
    Z_cls <- sweep(X_cls, 2L, mu_cls, FUN = "-")
    Z_class[idx_cls, ] <- sweep(Z_cls, 2L, sd_cls, FUN = "/")
  }
  
  class_cutoff <- 3
  idx <- which(is.finite(Z_class) & abs(Z_class) > class_cutoff, arr.ind = TRUE)
  row_ids <- idx[, "row"]
  col_ids <- idx[, "col"]
  flagged_z <- Z_class[idx]
  
  extreme_values <- data.frame(
    df[row_ids, meta_cols, drop = FALSE],
    metabolite = met_cols[col_ids],
    value_raw = X_raw[idx],
    value_log = X_log[idx],
    z_class = flagged_z,
    abs_z_class = abs(flagged_z),
    T2 = T2[row_ids],
    is_outlier_sample = is_outlier_sample[row_ids],
    z_threshold_used = rep(class_cutoff, nrow(idx)),
    stringsAsFactors = FALSE
  )
  
  extreme_values <- extreme_values[
    order(-extreme_values$abs_z_class, -extreme_values$T2),
    ,
    drop = FALSE
  ]
  rownames(extreme_values) <- NULL
  
  out_df <- df
  out_df$T2 <- T2
  out_df$is_outlier_sample <- is_outlier_sample
  out_df$used_in_fit <- used_in_fit
  
  z_class_df <- data.frame(
    df[, meta_cols, drop = FALSE],
    as.data.frame(Z_class, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  pca_plot <- NULL
  if (make_pca_plot && !is.null(pca)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' is required for PCA plotting but is not installed.")
    } else {
      plot_df <- data.frame(
        PC1 = scores_all[, 1L],
        PC2 = scores_all[, 2L],
        group = rep("Study Sample", length(idx_all))
      )
      plot_df$group[which(class_vec[idx_all] == qc_label)] <- "QC"
      plot_df$group[is_outlier_sample[idx_all]] <- "Potential Outlier"
      plot_df$group <- factor(
        plot_df$group,
        levels = c("QC", "Study Sample", "Potential Outlier")
      )
      
      eig <- eigen(cov_scores_reg, symmetric = TRUE)
      theta <- seq(0, 2 * pi, length.out = 300L)
      unit_circle <- rbind(cos(theta), sin(theta))
      transform <- eig$vectors %*% diag(sqrt(eig$values * cutoff), nrow = 2L)
      ellipse_pts <- sweep(transform %*% unit_circle, 1L, center_scores, FUN = "+")
      ellipse_df <- data.frame(PC1 = ellipse_pts[1L, ], PC2 = ellipse_pts[2L, ])
      
      pca_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2)) +
        ggplot2::geom_point(ggplot2::aes(color = group), alpha = 1, size = 4) +
        ggplot2::geom_path(
          data = ellipse_df,
          mapping = ggplot2::aes(x = PC1, y = PC2),
          inherit.aes = FALSE,
          linetype = "dashed"
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "QC" = "#F3C300",
            "Study Sample" = "#0067A5",
            "Potential Outlier" = "#E0003C"
          )
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 18, hjust = 0.5, face = "bold"),
          axis.title = ggplot2::element_text(size = 14, face = "bold"),
          axis.text = ggplot2::element_text(size = 14),
          legend.title = ggplot2::element_text(size = 14, face = "bold"),
          legend.text = ggplot2::element_text(size = 14)
        ) +
        ggplot2::labs(
          title = sprintf("PCA (PC1-PC2) with %g%% Mahalanobis ellipse", 100 * (1 - alpha)),
          #subtitle = "PCA fitted on non-QC samples. Ellipse estimated with QCs included.",
          color = "Group"
        )
    }
  }
  
  list(
    data = out_df,
    extreme_values = extreme_values,
    z_class = z_class_df,
    pca_plot = pca_plot,
    pc_loadings = pc_loadings,
    params = list(
      alpha = alpha,
      cutoff = cutoff,
      class_z_threshold = class_cutoff,
      log_transform = log_transform,
      log_offset = log_offset,
      p = ncol(X_scaled_all),
      class_col = class_col,
      qc_label = qc_label,
      min_complete = min_complete,
      drop_constant = drop_constant,
      const_tol = const_tol,
      ridge_factor = ridge_factor,
      retained_metabolites = met_cols
    )
  )
}