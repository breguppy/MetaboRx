#' @keywords internal
#' @noRd
.safe_loess_predict_x <- function(qc_x, qc_y, newx, span, degree) {
  ok <- is.finite(qc_x) & is.finite(qc_y) & qc_y > 0
  qx <- as.numeric(qc_x[ok])
  qy <- as.numeric(qc_y[ok])
  
  # require at least 2 points
  if (length(qx) < 2L) return(rep(1, length(newx)))
  
  # sort by x for approx + loess stability
  ord <- order(qx)
  qx <- qx[ord]
  qy <- qy[ord]
  
  n <- length(qx)
  
  # if too few QCs for stable loess, use interpolation
  # (this keeps behavior predictable for degree 0/1/2)
  if (n < (degree + 2L)) {
    return(stats::approx(qx, qy, xout = newx, rule = 2)$y)
  }
  
  # Respect the requested degree, but clamp to [0, 2]
  deg_req <- as.integer(degree)
  deg <- max(0L, min(2L, deg_req))
  
  # Span guardrail: ensure enough effective neighbors
  spn <- max(span, min(1, 8 / n))
  
  pred <- tryCatch({
    fit <- stats::loess(
      log(qy) ~ qx,
      span = spn,
      degree = deg,
      family = "symmetric",
      control = stats::loess.control(surface = "direct")
    )
    exp(stats::predict(fit, newdata = data.frame(qx = newx)))
  }, error = function(e) NA_real_)
  
  if (!is.numeric(pred) || all(!is.finite(pred))) {
    stats::approx(qx, qy, xout = newx, rule = 2)$y
  } else {
    pred
  }
}



loess_correction <- function(df, metab_cols, degree, span = 0.75, min_qc = 5) {
  df <- df[order(df$order), , drop = FALSE]
  if (!(identical(df$class[1], "QC") && identical(df$class[nrow(df)], "QC")))
    stop("First and last samples must be QCs.")
  
  qcid <- which(df$class == "QC")
  if (length(qcid) < min_qc) stop(sprintf("Need at least %d QC rows for LOESS.", min_qc))
  
  x_all <- suppressWarnings(as.numeric(df$order))
  if (any(!is.finite(x_all))) stop("order must be numeric and finite after sorting.")
  
  out <- df
  
  for (metab in metab_cols) {
    zero_mask <- is.finite(df[[metab]]) & df[[metab]] == 0
    qc_y <- df[[metab]][qcid]
    
    if (all(qc_y <= 0, na.rm = TRUE)) {
      out[[metab]] <- 0
      next
    }
    
    pred <- .safe_loess_predict_x(
      qc_x = x_all[qcid],
      qc_y = qc_y,
      newx = x_all,
      span = span,
      degree = degree
    )
    
    pred[!is.finite(pred) | pred <= 0] <- NA_real_
    corr <- as.numeric(df[[metab]]) / pred
    
    sf <- stats::median(corr[qcid], na.rm = TRUE)
    if (is.finite(sf) && sf > 0) corr <- corr / sf
    
    corr[!is.finite(corr) | corr < 0] <- NA_real_
    out[[metab]] <- corr
    out[[metab]][zero_mask] <- 0
  }
  
  # 99% NA fallback
  if (anyNA(out[metab_cols])) {
    for (metab in metab_cols) {
      x <- out[[metab]]
      if (mean(is.na(x)) >= 0.99) {
        mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
        x[is.na(x)] <- if (is.finite(mp)) mp else 0
        out[[metab]] <- x
      }
    }
    needs_knn <- anyNA(out[metab_cols]) && length(metab_cols) >= 2
    if (needs_knn) {
      kn <- impute::impute.knn(t(as.matrix(out[metab_cols])), rowmax = 1, colmax = 1, maxp = 15000)
      out[metab_cols] <- as.data.frame(t(kn$data))
    }
  }
  
  # final cleanup
  out[metab_cols] <- lapply(out[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    x[is.na(x)] <- if (is.finite(mp)) mp else 0
    x
  })
  
  out
}


#' Batch-wise LOESS correction (QC-based), using injection order as x
#'
#' Corrects each metabolite within each batch by dividing by a LOESS-smoothed
#' QC trend fit against `order` (not row index). Optionally normalizes the
#' corrected QC distribution globally per-metabolite so QC median is ~ 1.
#'
#' @param df        Data frame containing metadata columns: batch, class, order
#'                 plus metabolite columns.
#' @param metab_cols Character vector of metabolite column names.
#' @param span      LOESS span parameter.
#' @param degree    LOESS degree (max 2). Internally reduced for small QC counts.
#' @param min_qc    Minimum QC rows required within each batch to attempt LOESS.
#'
#' @return Data frame with corrected metabolite columns.
#'
#' @keywords internal
#' @noRd
bw_loess_correction <- function(df, metab_cols, span = 0.75, degree = 2, min_qc = 5) {
  if (!all(c("batch", "class", "order") %in% names(df))) {
    stop("df must contain columns: batch, class, order.")
  }
  if (!all(metab_cols %in% names(df))) {
    missing <- setdiff(metab_cols, names(df))
    stop(sprintf("Missing metabolite columns: %s", paste(missing, collapse = ", ")))
  }
  
  out <- df
  
  for (metab in metab_cols) {
    # preserve original exact zeros (common for missing/below-LOD encoding)
    zero_mask <- is.finite(df[[metab]]) & df[[metab]] == 0
    
    for (b in unique(df$batch)) {
      b_idx <- which(df$batch == b)
      sub   <- df[b_idx, , drop = FALSE]
      
      if (!(identical(sub$class[1], "QC") && identical(sub$class[nrow(sub)], "QC"))) {
        stop(sprintf("Batch '%s' must start and end with QC.", b))
      }
      
      qcid <- which(sub$class == "QC")
      if (length(qcid) < min_qc) {
        warning(sprintf(
          "Skipping batch '%s' for '%s': only %d QC rows (< %d).",
          b, metab, length(qcid), min_qc
        ))
        next
      }
      
      qc_y <- sub[[metab]][qcid]
      if (all(qc_y <= 0, na.rm = TRUE)) {
        out[[metab]][b_idx] <- 0
        next
      }
      
      x_sub <- suppressWarnings(as.numeric(sub$order))
      if (any(!is.finite(x_sub))) {
        stop(sprintf("Non-numeric or non-finite `order` detected in batch '%s'.", b))
      }
      
      pred <- .safe_loess_predict_x(
        qc_x   = x_sub[qcid],
        qc_y   = qc_y,
        newx   = x_sub,
        span   = span,
        degree = degree
      )
      
      pred[!is.finite(pred) | pred <= 0] <- NA_real_
      corr <- as.numeric(sub[[metab]]) / pred
      corr[!is.finite(corr) | corr < 0] <- NA_real_
      
      out[[metab]][b_idx] <- corr
    }
    
    # Global re-anchoring per metabolite: QC median -> 1
    qc_all <- out$class == "QC" & is.finite(out[[metab]]) & out[[metab]] > 0
    gsf <- suppressWarnings(stats::median(out[[metab]][qc_all], na.rm = TRUE))
    if (is.finite(gsf) && gsf > 0) {
      out[[metab]] <- out[[metab]] / gsf
    }
    
    # preserve exact zeros
    out[[metab]][zero_mask] <- 0
  }
  
  # If there are NAs, do:
  # 1) per-metabolite "almost all NA" fallback
  # 2) kNN impute only if >= 2 metabolites (otherwise it can invent values)
  if (anyNA(out[metab_cols])) {
    for (metab in metab_cols) {
      x <- out[[metab]]
      if (mean(is.na(x)) >= 0.99) {
        mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
        x[is.na(x)] <- if (is.finite(mp)) mp else 0
        out[[metab]] <- x
      }
    }
    
    needs_knn <- anyNA(out[metab_cols]) && length(metab_cols) >= 2
    if (needs_knn) {
      kn <- impute::impute.knn(
        t(as.matrix(out[metab_cols])),
        rowmax = 1,
        colmax = 1,
        maxp   = 15000
      )
      out[metab_cols] <- as.data.frame(t(kn$data))
    }
  }
  
  # Final cleanup: enforce non-negative finite values, fill remaining NA with
  # smallest positive, else 0.
  out[metab_cols] <- lapply(out[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    x[is.na(x)] <- if (is.finite(mp)) mp else 0
    x
  })
  
  out
}
