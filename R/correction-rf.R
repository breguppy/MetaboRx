#' @keywords internal
#' @noRd
.safe_baseline_predict <- function(qcid, qc_y, newx) {
  ok <- is.finite(qc_y)
  # flat baseline
  if (sum(ok) < 2L) {
    return(rep(1, length(newx)))
  }
  # Otherwise linear interpolation
  qx <- qcid[ok]
  qy <- qc_y[ok]
  stats::approx(qx, qy, xout = newx, rule = 2)$y
}

#' @keywords internal
#' @noRd
.safe_rf_predict <- function(qc_order, qc_y, new_order, ntree, quiet = TRUE) {
  ok <- is.finite(qc_order) & is.finite(qc_y)
  if (sum(ok) < 2L || length(unique(qc_y[ok])) < 2L) {
    return(.safe_baseline_predict(qc_order, qc_y, new_order))
  }

  rf_call <- function() {
    fit <- randomForest::randomForest(
      x = data.frame(order = qc_order[ok]),
      y = qc_y[ok],
      ntree = ntree
    )
    stats::predict(fit, newdata = data.frame(order = new_order))
  }

  pred <- tryCatch(
    {
      if (quiet) {
        withCallingHandlers(
          rf_call(),
          warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("five or fewer unique values", msg, ignore.case = TRUE)) {
              invokeRestart("muffleWarning")
            }
          }
        )
      } else {
        rf_call()
      }
    },
    error = function(e) NA_real_
  )

  if (!is.numeric(pred) || all(!is.finite(pred))) {
    .safe_baseline_predict(qc_order, qc_y, new_order)
  } else {
    pred
  }
}

#' Random forest correction
#'
#' @param df Data frame containing class, order, and metabolite columns.
#' @param metab_cols Character vector of metabolite column names.
#' @param ntree Number of trees.
#' @param seed Optional random seed.
#' @param min_qc Minimum QC count for fitting a random forest.
#'   Below this threshold, use baseline interpolation.
#' @param met_scaling Character. If "original", rescale each corrected
#'   metabolite so its QC median matches its raw QC median.
#'   Only finite, positive raw QC values with valid corrected values
#'   are used to calculate the scaling factor.
#'
#' @keywords internal
#' @noRd
rf_correction <- function(
    df,
    metab_cols,
    ntree = 500,
    seed = NULL,
    min_qc = 5,
    met_scaling = "one"
) {
  if (!is.null(seed)) set.seed(seed)
  
  df <- df[order(df$order), , drop = FALSE]
  
  if (!(identical(df$class[1], "QC") &&
        identical(df$class[nrow(df)], "QC"))) {
    stop("First and last samples must be QCs.")
  }
  
  out <- df
  qcid <- which(df$class == "QC")
  newx <- df$order
  
  use_baseline <- length(qcid) < min_qc
  if (use_baseline) {
    warning(
      sprintf(
        "Fewer than %d QC samples; using baseline interpolation.",
        min_qc
      )
    )
  }
  
  for (metab in metab_cols) {
    raw_y <- as.numeric(df[[metab]])
    zero_mask <- is.finite(raw_y) & raw_y == 0
    qc_y <- raw_y[qcid]
    
    if (all(qc_y <= 0, na.rm = TRUE)) {
      out[[metab]] <- 0
      next
    }
    
    pred <- if (use_baseline) {
      .safe_baseline_predict(
        qcid = newx[qcid],
        qc_y = qc_y,
        newx = newx
      )
    } else {
      .safe_rf_predict(
        qc_order = newx[qcid],
        qc_y = qc_y,
        new_order = newx,
        ntree = ntree
      )
    }
    
    pred[!is.finite(pred) | pred <= 0] <- NA_real_
    
    corr <- raw_y / pred
    corr[!is.finite(corr) | corr < 0] <- NA_real_
    
    valid_qc <- is.finite(corr[qcid]) & corr[qcid] > 0
    
    if (any(valid_qc)) {
      corr <- corr / stats::median(corr[qcid][valid_qc])
    }
    
    if (met_scaling == "original") {
      corrected_qc <- corr[qcid]
      
      # Use the same valid QC rows for both medians.
      valid_qc <- is.finite(qc_y) &
        qc_y > 0 &
        is.finite(corrected_qc) &
        corrected_qc > 0
      
      if (any(valid_qc)) {
        raw_qc_median <- stats::median(qc_y[valid_qc])
        corrected_qc_median <- stats::median(
          corrected_qc[valid_qc]
        )
        
        scale_factor <- raw_qc_median / corrected_qc_median
        
        if (is.finite(scale_factor) && scale_factor > 0) {
          corr <- corr * scale_factor
        } else {
          warning(
            sprintf("Could not preserve QC scale for '%s'.", metab)
          )
        }
      } else {
        warning(
          sprintf(
            "Could not preserve QC scale for '%s': no valid positive QC pairs.",
            metab
          )
        )
      }
    }
    
    corr[!is.finite(corr) | corr < 0] <- NA_real_
    corr[zero_mask] <- 0
    out[[metab]] <- corr
  }
  
  .cleanup_corrected_metabolites(out, metab_cols)
}

#' Batch-wise random forest correction
#'
#' @keywords internal
#' @noRd
bw_rf_correction <- function(df, metab_cols, ntree = 500, seed = NULL, min_qc = 5) {
  if (!is.null(seed)) set.seed(seed)
  out <- df
  batches <- unique(df$batch)
  rows_by_batch <- lapply(batches, function(b) {
    idx_all <- which(df$batch == b)
    idx_all[order(df$order[idx_all])]
  })

  for (metab in metab_cols) {
    zero_mask <- is.finite(df[[metab]]) & df[[metab]] == 0
    for (i in seq_along(batches)) {
      b <- batches[[i]]
      b_idx <- rows_by_batch[[i]]
      sub <- df[b_idx, , drop = FALSE]

      # require QC at start and end within batch
      if (!(identical(sub$class[1], "QC") && identical(sub$class[nrow(sub)], "QC"))) {
        stop(sprintf("Batch '%s' must start and end with QC.", b))
      }

      qcid <- which(sub$class == "QC")
      qcy <- sub[[metab]][qcid]
      if (all(qcy <= 0, na.rm = TRUE)) {
        out[[metab]][b_idx] <- 0
        next
      }
      if (length(qcid) < min_qc) {
        warning(sprintf("Skipping %s in batch %s - too few QC samples; using baseline interpolation.", metab, b))
      }

      pred <- .safe_rf_predict(sub$order[qcid], sub[[metab]][qcid],
        new_order = sub$order, ntree = ntree
      )
      pred[!is.finite(pred) | pred <= 0] <- NA_real_
      corr <- as.numeric(sub[[metab]]) / pred
      corr[!is.finite(corr) | corr < 0] <- NA_real_
      out[[metab]][b_idx] <- corr
      out[[metab]][zero_mask] <- 0
    }
  }

  .cleanup_corrected_metabolites(out, metab_cols)
}
