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
#' @keywords internal
#' @noRd
rf_correction <- function(df, metab_cols, ntree = 500, seed = NULL, min_qc = 5) {
  if (!is.null(seed)) set.seed(seed)

  # sort and require QC at start and end
  df <- df[order(df$order), , drop = FALSE]
  if (!(identical(df$class[1], "QC") && identical(df$class[nrow(df)], "QC"))) {
    stop("First and last samples must be QCs.")
  }

  out <- df
  qcid <- which(df$class == "QC")
  if (length(qcid) < 2L) warning("Very few QC rows; falling back to baseline where needed.")
  newx <- df$order

  for (metab in metab_cols) {
    zero_mask <- is.finite(df[[metab]]) & df[[metab]] == 0
    qc_y <- df[[metab]][qcid]
    if (all(qc_y <= 0, na.rm = TRUE)) {
      out[[metab]] <- 0
      next
    }
    if (length(qcid) < min_qc) {
      warning(paste("Skipping", metab, "- too few QC samples for RF; using baseline interpolation."))
    }
    pred <- .safe_rf_predict(df$order[qcid], qc_y, newx, ntree = ntree)
    pred[!is.finite(pred) | pred <= 0] <- NA_real_
    corr <- as.numeric(df[[metab]]) / pred
    corr[!is.finite(corr) | corr < 0] <- NA_real_
    out[[metab]] <- corr
    out[[metab]][zero_mask] <- 0
  }

  # final cleanup: smallest positive per metabolite, else 0
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
