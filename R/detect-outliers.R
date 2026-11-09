#' Detect sample and metabolite outliers
#'
#' Fits PCA on scaled, complete non-QC rows and estimates the PC1-PC2
#' Mahalanobis reference from all complete projected rows, including QCs.
#'
#' Groups with fewer than 8 values were tested with a two-sided Grubbs-type
#' test (alpha = 0.005) that used the within-group SD pooled across all groups
#' for that metabolite (externally studentized one-way ANOVA residuals,
#' Bonferroni-adjusted within each group); larger groups were tested with
#' ROUT (Q = 0.1% by default). Both were applied to log-transformed data,
#' and flagged values were not excluded.
#' Selection is independent of the sample's Mahalanobis flag.
#'
#' Missing values are imputed for PCA. Supply observed_df containing corrected
#' intensities with original missing cells restored to NA when df has already
#' been imputed. Without observed_df, only missing cells still in df can be
#' excluded. Output rows preserve the input order.
#'
#' @param observed_df Optional corrected data frame with observed intensities
#'   only. Must have the same unique row keys, class labels and metabolites.
#' @param row_key_cols Columns uniquely identifying samples in both inputs.
#'
#' Legacy threshold arguments are accepted for compatibility and ignored.
#' Legacy Z-score fields contain signed or absolute test statistics.
#' @param rout_q Target false discovery rate per metabolite and class.
#' `used_in_fit` identifies rows used to fit PCA.
#'
#' @keywords internal
#' @noRd
detect_outliers <- function(
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
    make_pca_plot = TRUE,
    observed_df = NULL,
    row_key_cols = c("sample", "batch"),
    rout_q = 0.001
) {
  if (!is.numeric(rout_q) || length(rout_q) != 1L ||
      !is.finite(rout_q) || rout_q <= 0 || rout_q >= 1) {
    stop("rout_q must be a single number strictly between 0 and 1.")
  }
  
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
  
  if (length(row_key_cols) == 0L || anyDuplicated(row_key_cols)) {
    stop("row_key_cols must contain distinct identifying column names.")
  }
  
  # Length prefixes make keys unambiguous even when labels contain separators.
  row_keys <- function(x) {
    missing_keys <- setdiff(row_key_cols, names(x))
    if (length(missing_keys) > 0L) {
      stop("Missing row key columns: ", paste(missing_keys, collapse = ", "))
    }
    if (anyNA(x[, row_key_cols, drop = FALSE])) {
      stop("Row key columns must not contain missing values.")
    }
    parts <- lapply(x[row_key_cols], function(column) {
      values <- as.character(column)
      paste0(nchar(values, type = "bytes"), ":", values)
    })
    keys <- do.call(paste0, parts)
    if (anyDuplicated(keys)) {
      stop("row_key_cols must uniquely identify every row.")
    }
    keys
  }
  
  input_keys <- row_keys(df)
  align_rows <- function(x) {
    keys <- row_keys(x)
    if (length(keys) != length(input_keys) || !setequal(keys, input_keys)) {
      stop("Input and reference data must contain the same row keys.")
    }
    x <- x[match(input_keys, keys), , drop = FALSE]
    if (!class_col %in% names(x) ||
        !identical(as.character(x[[class_col]]), as.character(df[[class_col]]))) {
      stop("Class labels must agree for matching row keys.")
    }
    x
  }
  
  observed_source <- if (is.null(observed_df)) df else align_rows(observed_df)
  if (!all(met_cols %in% names(observed_source)) ||
      !all(vapply(observed_source[met_cols], is.numeric, logical(1)))) {
    stop("observed_df must contain all numeric metabolite columns from df.")
  }
  X_observed <- as.matrix(observed_source[, met_cols, drop = FALSE])
  # Cells still missing in the PCA input cannot be treated as observed.
  X_observed[!is.finite(as.matrix(df[, met_cols, drop = FALSE]))] <- NA_real_
  
  if (anyNA(df[, met_cols, drop = FALSE])) {
    results <- impute_missing(df, met_cols, p$qcImputeM, p$samImputeM)
    # Restore input order before combining PCA data, observed values and metadata.
    df <- align_rows(results$df)
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
  
  # Value screening always uses log2 intensities, independently of PCA.
  X_observed_log <- log2(X_observed + log_offset)
  X_observed_log[!is.finite(X_observed_log)] <- NA_real_
  
  value_scores <- matrix(
    NA_real_,
    nrow = nrow(df),
    ncol = ncol(X_observed_log),
    dimnames = dimnames(X_observed_log)
  )
  adjusted_p <- rout_scores <- value_scores
  value_method <- matrix(NA_character_, nrow(df), ncol(X_observed_log))
  grubbs_alpha <- 0.005
  min_rout_observations <- 8L
  scoring_groups <- split(seq_along(class_vec), class_vec, drop = TRUE)
  
  for (col in seq_len(ncol(X_observed_log))) {
    column <- X_observed_log[, col]
    observed_groups <- lapply(scoring_groups, function(rows) {
      rows[is.finite(column[rows])]
    })
    counts <- lengths(observed_groups)
    pooled_df <- sum(counts[counts > 0L] - 1L)
    group_sse <- vapply(observed_groups, function(rows) {
      values <- column[rows]
      sum((values - mean(values))^2)
    }, numeric(1))
    
    for (group in seq_along(observed_groups)) {
      observed_rows <- observed_groups[[group]]
      count <- counts[group]
      values <- column[observed_rows]
      
      if (count < min_rout_observations) {
        if (count < 2L || pooled_df <= 1L) {
          next
        }
        # Recompute the small group's SSE after deletion to avoid cancellation.
        other_sse <- sum(group_sse[-group])
        deleted_sse <- vapply(seq_len(count), function(i) {
          remaining <- values[-i]
          other_sse + sum((remaining - mean(remaining))^2)
        }, numeric(1))
        deleted_df <- pooled_df - 1L
        denominator <- sqrt(deleted_sse / deleted_df * (1 - 1 / count))
        scores <- rep(NA_real_, count)
        valid <- is.finite(denominator) & denominator > 0
        scores[valid] <- (values[valid] - mean(values)) / denominator[valid]
        p_values <- 2 * stats::pt(-abs(scores), df = deleted_df)
        value_scores[observed_rows, col] <- scores
        adjusted_p[observed_rows, col] <- stats::p.adjust(
          p_values, method = "bonferroni", n = count
        )
        value_method[observed_rows, col] <- "Grubbs-type"
        next
      }
      
      # Intercept-only Lorentzian fit with adaptive residual scale.
      values <- values - stats::median(values)
      span <- max(abs(values))
      if (!is.finite(span) || span == 0) {
        next
      }
      values <- values / span
      
      center <- 0
      converged <- FALSE
      for (iteration in seq_len(200L)) {
        residuals <- values - center
        rsdr <- stats::quantile(
          abs(residuals), 0.6827, names = FALSE, type = 7
        ) * count / (count - 1L)
        if (!is.finite(rsdr) || rsdr <= 0) {
          break
        }
        
        weights <- 1 / (1 + (residuals / rsdr)^2)
        next_center <- stats::weighted.mean(values, weights)
        if (abs(next_center - center) <= sqrt(.Machine$double.eps) * rsdr) {
          converged <- TRUE
          break
        }
        center <- next_center
      }
      if (!converged) {
        next
      }
      
      scores <- residuals / rsdr
      rout_scores[observed_rows, col] <- scores
      value_scores[observed_rows, col] <- scores
      value_method[observed_rows, col] <- "ROUT"
      # Round the 30% tail up, retaining all observations in the BH denominator.
      tested <- head(order(abs(scores), decreasing = TRUE), ceiling(0.30 * count))
      p_values <- 2 * stats::pt(-abs(scores[tested]), df = count - 1L)
      adjusted_p[observed_rows[tested], col] <- stats::p.adjust(
        p_values, method = "BH", n = count
      )
    }
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
  
  value_cutoff <- ifelse(value_method == "ROUT", rout_q, grubbs_alpha)
  extreme_mask <- !is.na(adjusted_p) & adjusted_p < value_cutoff
  idx <- which(extreme_mask, arr.ind = TRUE)
  row_ids <- idx[, "row"]
  col_ids <- idx[, "col"]
  flagged_z <- value_scores[idx]
  
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
    method = value_method[idx],
    test_statistic = flagged_z,
    adjusted_p = adjusted_p[idx],
    p_cutoff = value_cutoff[idx],
    rout_score = rout_scores[idx],
    rout_adjusted_p = replace(adjusted_p[idx], value_method[idx] != "ROUT", NA_real_),
    z_threshold_used = rep(NA_real_, nrow(idx)),
    flag_fraction = rep(NA_real_, nrow(idx)),
    stringsAsFactors = FALSE
  )
  
  extreme_values <- extreme_values[
    order(extreme_values$adjusted_p, -abs(extreme_values$test_statistic),
          -extreme_values$T2),
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
    as.data.frame(value_scores, check.names = FALSE),
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
            "Study Sample" = "#6180D9",  #"#2A4EB3", # "#4169E1",
            "Potential Outlier" = "#FF0328"# "#D2232A"  #"#E0003C"
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
    rout_scores = data.frame(
      df[, meta_cols, drop = FALSE],
      as.data.frame(rout_scores, check.names = FALSE),
      check.names = FALSE
    ),
    value_scores = z_global_df,
    z_global = z_global_df,
    z_class = z_global_df,
    pca_plot = pca_plot,
    pc_loadings = pc_loadings,
    params = list(
      observed_reference_supplied = !is.null(observed_df),
      row_key_cols = row_key_cols,
      alpha = alpha,
      cutoff = cutoff,
      value_z_threshold = NA_real_,
      rout_q = rout_q,
      flag_fraction = NA_real_,
      min_score_observations = 2L,
      min_rout_observations = min_rout_observations,
      grubbs_alpha = grubbs_alpha,
      flagging_rule = "Grubbs-type for fewer than 8 observed values; ROUT otherwise",
      z_score_reference = "Class-specific centers including QC; Grubbs-type scale pooled across classes",
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