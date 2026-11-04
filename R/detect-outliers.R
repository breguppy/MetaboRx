#' Detect sample and metabolite outliers
#'
#' Fits PCA on scaled, complete non-QC rows and estimates the PC1-PC2
#' Mahalanobis reference from all complete projected rows, including QCs.
#'
#' Computes metabolite Z-scores separately for pooled study samples and QCs
#' using finite observed values on the selected transformation scale.
#' Requires an absolute Z-score greater than 2, then reports the largest
#' 10% of eligible values across metabolites within each pool, rounding
#' counts down and breaking ties by row then column.
#' Selection is independent of the sample's Mahalanobis flag.
#'
#' Missing values are imputed for PCA. Imputed values are excluded from
#' the metabolite Z-score calculation.
#'
#' Legacy threshold arguments are accepted for compatibility and ignored.
#' `z_class` and `abs_z_class` are compatibility aliases for pooled scores.
#' `used_in_fit` identifies rows used to fit PCA.
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
    stop(
      "Missing metadata columns in df: ",
      paste(missing_meta, collapse = ", ")
    )
  }
  if (!class_col %in% names(df)) {
    stop("class_col '", class_col, "' not found in df.")
  }
  
  candidate_cols <- setdiff(names(df), c(meta_cols, class_col))
  met_cols <- candidate_cols[
    vapply(df[candidate_cols], is.numeric, logical(1))
  ]
  if (length(met_cols) == 0L) {
    stop("No numeric metabolite columns found.")
  }
  
  # Preserve observed values for the metabolite screen.
  X_observed <- as.matrix(df[, met_cols, drop = FALSE])
  
  if (anyNA(X_observed)) {
    results <- impute_missing(df, met_cols, p$qcImputeM, p$samImputeM)
    df <- results$df
  }
  
  X_raw <- as.matrix(df[, met_cols, drop = FALSE])
  X_log <- if (log_transform) log2(X_raw + log_offset) else X_raw
  
  class_vec <- df[[class_col]]
  nonqc_mask <- !is.na(class_vec) & class_vec != qc_label
  complete_nonqc <- nonqc_mask &
    rowSums(!is.finite(X_log)) == 0L
  X_fit <- X_log[complete_nonqc, , drop = FALSE]
  
  if (sum(complete_nonqc) < min_complete) {
    stop("Too few complete non-QC rows to estimate PCA.")
  }
  
  # Pooled metabolite Z-scores include every numeric metabolite,
  # including metabolites subsequently excluded from PCA.
  X_observed_log <- if (log_transform) {
    log2(X_observed + log_offset)
  } else {
    X_observed
  }
  X_observed_log[!is.finite(X_observed_log)] <- NA_real_
  
  # Compute separate pooled Z-scores for study samples and QCs.
  # Missing class labels remain unscored.
  Z_global <- matrix(
    NA_real_,
    nrow = nrow(df),
    ncol = ncol(X_observed_log),
    dimnames = dimnames(X_observed_log)
  )
  
  scoring_groups <- list(
    study = which(!is.na(class_vec) & class_vec != qc_label),
    qc = which(!is.na(class_vec) & class_vec == qc_label)
  )
  
  for (rows in scoring_groups) {
    if (length(rows) == 0L) {
      next
    }
    
    X_group <- X_observed_log[rows, , drop = FALSE]
    group_mean <- colMeans(X_group, na.rm = TRUE)
    group_sd <- apply(X_group, 2L, stats::sd, na.rm = TRUE)
    group_sd[!is.finite(group_sd) | group_sd == 0] <- NA_real_
    
    Z_group <- sweep(X_group, 2L, group_mean, FUN = "-")
    Z_global[rows, ] <- sweep(Z_group, 2L, group_sd, FUN = "/")
  }
  
  # Retain all metabolite names for the value-level screen.
  screened_metabolites <- met_cols
  
  if (drop_constant) {
    variances <- apply(X_fit, 2L, stats::var, na.rm = TRUE)
    const_mask <- !is.finite(variances) | variances <= const_tol
    
    if (any(const_mask)) {
      message(
        "Excluding ", sum(const_mask),
        " metabolite(s) from PCA with near-zero variance among non-QC: ",
        paste(met_cols[const_mask], collapse = ", ")
      )
      
      keep_mask <- !const_mask
      met_cols <- met_cols[keep_mask]
      X_log <- X_log[, keep_mask, drop = FALSE]
      X_fit <- X_fit[, keep_mask, drop = FALSE]
    }
  }
  
  n <- nrow(df)
  T2 <- rep(NA_real_, n)
  used_in_fit <- rep(FALSE, n)
  cutoff <- stats::qchisq(1 - alpha, df = 2L)
  
  pca <- NULL
  pc_loadings <- NULL
  pca_plot <- NULL
  
  if (length(met_cols) >= 2L) {
    X_fit_scaled <- scale(X_fit)
    center_scaled <- attr(X_fit_scaled, "scaled:center")
    scale_scaled <- attr(X_fit_scaled, "scaled:scale")
    scale_scaled[!is.finite(scale_scaled) | scale_scaled == 0] <- NA_real_
    
    X_scaled_all <- sweep(X_log, 2L, center_scaled, FUN = "-")
    X_scaled_all <- sweep(X_scaled_all, 2L, scale_scaled, FUN = "/")
    
    complete_scaled <- rowSums(!is.finite(X_scaled_all)) == 0L
    idx_all <- which(complete_scaled)
    idx_fit <- which(nonqc_mask & complete_scaled)
    
    if (length(idx_fit) >= 3L) {
      pca <- stats::prcomp(
        X_scaled_all[idx_fit, , drop = FALSE],
        center = FALSE,
        scale. = FALSE
      )
      
      rotation <- pca$rotation[, 1:2, drop = FALSE]
      scores_all <- X_scaled_all[idx_all, , drop = FALSE] %*% rotation
      
      center_scores <- colMeans(scores_all)
      cov_scores <- stats::cov(scores_all)
      mean_diag_pc <- mean(diag(cov_scores))
      
      if (!is.finite(mean_diag_pc) || mean_diag_pc <= 0) {
        mean_diag_pc <- 1
      }
      
      cov_scores_reg <- cov_scores +
        diag(ridge_factor * mean_diag_pc, 2L)
      
      T2[idx_all] <- stats::mahalanobis(
        scores_all,
        center_scores,
        cov_scores_reg
      )
      used_in_fit[idx_fit] <- TRUE
      
      pc_loadings <- data.frame(
        metabolite = met_cols,
        PC1 = rotation[, 1L],
        PC2 = rotation[, 2L],
        stringsAsFactors = FALSE
      )
    }
  }
  
  is_outlier_sample <- !is.na(T2) & T2 > cutoff
  
  value_z_threshold <- 2
  flag_fraction <- 0.10
  extreme_mask <- matrix(
    FALSE,
    nrow = nrow(Z_global),
    ncol = ncol(Z_global),
    dimnames = dimnames(Z_global)
  )
  
  for (rows in scoring_groups) {
    if (length(rows) == 0L) {
      next
    }
    
    group_z <- Z_global[rows, , drop = FALSE]
    usable <- which(
      is.finite(group_z) & abs(group_z) > value_z_threshold,
      arr.ind = TRUE
    )
    n_flag <- floor(nrow(usable) * flag_fraction)
    
    if (n_flag == 0L) {
      next
    }
    
    ranked <- order(
      -abs(group_z[usable]),
      usable[, "row"],
      usable[, "col"]
    )
    selected <- usable[ranked[seq_len(n_flag)], , drop = FALSE]
    extreme_mask[cbind(rows[selected[, "row"]], selected[, "col"])] <- TRUE
  }
  
  idx <- which(extreme_mask, arr.ind = TRUE)
  row_ids <- idx[, "row"]
  col_ids <- idx[, "col"]
  flagged_z <- Z_global[idx]
  
  extreme_values <- data.frame(
    df[row_ids, meta_cols, drop = FALSE],
    metabolite = screened_metabolites[col_ids],
    value_raw = X_observed[idx],
    value_log = X_observed_log[idx],
    z_global = flagged_z,
    abs_z_global = abs(flagged_z),
    z_class = flagged_z,
    abs_z_class = abs(flagged_z),
    T2 = T2[row_ids],
    is_outlier_sample = is_outlier_sample[row_ids],
    z_threshold_used = rep(value_z_threshold, nrow(idx)),
    flag_fraction = rep(flag_fraction, nrow(idx)),
    stringsAsFactors = FALSE
  )
  
  extreme_values <- extreme_values[
    order(-extreme_values$abs_z_global, -extreme_values$T2),
    ,
    drop = FALSE
  ]
  rownames(extreme_values) <- NULL
  
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
  
  if (make_pca_plot && !is.null(pca)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning(
        "Package 'ggplot2' is required for PCA plotting but is not installed."
      )
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
      transform <- eig$vectors %*%
        diag(sqrt(eig$values * cutoff), nrow = 2L)
      ellipse_pts <- sweep(
        transform %*% unit_circle,
        1L,
        center_scores,
        FUN = "+"
      )
      ellipse_df <- data.frame(
        PC1 = ellipse_pts[1L, ],
        PC2 = ellipse_pts[2L, ]
      )
      
      pca_plot <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = PC1, y = PC2)
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = group),
          alpha = 1,
          size = 4
        ) +
        ggplot2::geom_path(
          data = ellipse_df,
          mapping = ggplot2::aes(x = PC1, y = PC2),
          inherit.aes = FALSE,
          linetype = "dashed"
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "QC" = "#F3C300",
            "Study Sample" =  "#2A4EB3", # "#4169E1",
            "Potential Outlier" = "#FF0300"# "#D2232A"  #"#E0003C"
          )
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 18, hjust = 0.5, face = "bold"
          ),
          axis.title = ggplot2::element_text(size = 14, face = "bold"),
          axis.text = ggplot2::element_text(size = 14),
          legend.title = ggplot2::element_text(size = 14, face = "bold"),
          legend.text = ggplot2::element_text(size = 14)
        ) +
        ggplot2::labs(
          title = sprintf(
            "PCA (PC1-PC2) with %g%% Mahalanobis ellipse",
            100 * (1 - alpha)
          ),
          color = "Group"
        )
    }
  }
  
  list(
    data = out_df,
    extreme_values = extreme_values,
    z_global = z_global_df,
    z_class = z_global_df,
    pca_plot = pca_plot,
    pc_loadings = pc_loadings,
    params = list(
      alpha = alpha,
      cutoff = cutoff,
      value_z_threshold = value_z_threshold,
      flag_fraction = flag_fraction,
      flagging_rule = "Largest 10% of values with absolute Z-score greater than 2 within each scoring pool",
      z_score_reference = "Non-QC samples pooled together; QCs scored separately",
      log_transform = log_transform,
      log_offset = log_offset,
      p = length(met_cols),
      class_col = class_col,
      qc_label = qc_label,
      min_complete = min_complete,
      drop_constant = drop_constant,
      const_tol = const_tol,
      ridge_factor = ridge_factor,
      retained_metabolites = met_cols,
      screened_metabolites = screened_metabolites
    )
  )
}