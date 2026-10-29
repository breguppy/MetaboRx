#' Safe adaptive Nadaraya-Watson prediction
#'
#' @param qc_x numeric vector
#' @param qc_y numeric vector
#' @param newx numeric vector
#' @param span numeric scalar in (0, 1]
#' @param kernel character scalar
#'
#' @return numeric vector
#' @keywords internal
#' @noRd
.safe_nw_predict_x <- function(qc_x, qc_y, newx, span = 0.75, kernel = "gaussian") {
  ok <- is.finite(qc_x) & is.finite(qc_y) & qc_y > 0
  qx <- as.numeric(qc_x[ok])
  qy <- as.numeric(qc_y[ok])

  if (length(qx) < 2L) {
    return(rep(1, length(newx)))
  }

  ord <- order(qx)
  qx <- qx[ord]
  qy <- qy[ord]

  n <- length(qx)
  k <- max(2L, min(n, ceiling(span * n)))

  kernel_fun <- switch(kernel,
    gaussian = function(u) stats::dnorm(u),
    tricube = function(u) {
      out <- (1 - abs(u)^3)^3
      out[abs(u) >= 1] <- 0
      out
    },
    stop("Unsupported kernel: ", kernel)
  )

  pred <- vapply(
    X = newx,
    FUN.VALUE = numeric(1),
    FUN = function(x0) {
      d <- abs(qx - x0)
      h <- sort(d, partial = k)[k]

      if (!is.finite(h) || h <= 0) {
        h <- max(d[d > 0], na.rm = TRUE)
      }
      if (!is.finite(h) || h <= 0) {
        return(stats::median(qy, na.rm = TRUE))
      }

      u <- (x0 - qx) / h
      w <- kernel_fun(u)
      sw <- sum(w, na.rm = TRUE)

      if (!is.finite(sw) || sw <= 0) {
        return(stats::approx(qx, qy, xout = x0, rule = 2)$y)
      }

      sum(w * qy, na.rm = TRUE) / sw
    }
  )

  bad <- !is.finite(pred) | pred <= 0
  if (any(bad)) {
    pred_fallback <- stats::approx(qx, qy, xout = newx, rule = 2)$y
    pred[bad] <- pred_fallback[bad]
  }

  pred
}

#' Nadaraya-Watson QC-based correction
#'
#' @param df Data frame containing class, order, and metabolite columns.
#' @param metab_cols Character vector of metabolite column names.
#' @param met_scaling Either "one" or "original". Final scaling uses
#'   matched finite, positive raw and corrected QC pairs identified
#'   before imputation. "one" sets their corrected median to 1;
#'   "original" sets it to their raw median.
#' @param span Relative smoothing span in (0, 1].
#' @param min_qc Minimum number of QC rows required.
#' @param kernel Either "gaussian" or "tricube".
#'
#' @return Corrected data frame.
#' @keywords internal
#' @noRd
nw_correction <- function(
    df,
    metab_cols,
    met_scaling = "one",
    span = 0.75,
    min_qc = 3,
    kernel = "gaussian"
) {
  met_scaling <- match.arg(met_scaling, c("one", "original"))
  kernel <- match.arg(kernel, c("gaussian", "tricube"))
  
  if (nrow(df) == 0L) {
    stop("df must contain at least one row.")
  }
  
  x_all <- suppressWarnings(as.numeric(df$order))
  
  if (any(!is.finite(x_all))) {
    stop("order must be numeric and finite.")
  }
  
  ord <- order(x_all)
  df <- df[ord, , drop = FALSE]
  x_all <- x_all[ord]
  
  if (!(isTRUE(df$class[1] == "QC") &&
        isTRUE(df$class[nrow(df)] == "QC"))) {
    stop("First and last samples must be QCs.")
  }
  
  qcid <- which(df$class == "QC")
  
  if (length(qcid) < min_qc) {
    stop(
      sprintf(
        "Need at least %d QC rows for local constant correction.",
        min_qc
      )
    )
  }
  
  out <- df
  
  for (metab in metab_cols) {
    raw_y <- as.numeric(df[[metab]])
    qc_y <- raw_y[qcid]
    
    if (!any(is.finite(qc_y) & qc_y > 0)) {
      out[[metab]] <- 0
      next
    }
    
    pred <- .safe_nw_predict_x(
      qc_x = x_all[qcid],
      qc_y = qc_y,
      newx = x_all,
      span = span,
      kernel = kernel
    )
    
    pred[!is.finite(pred) | pred <= 0] <- NA_real_
    
    corr <- raw_y / pred
    corr[!is.finite(corr) | corr < 0] <- NA_real_
    
    out[[metab]] <- corr
  }
  
  .finalize_qc_correction(
    out = out,
    raw_df = df,
    metab_cols = metab_cols,
    qcid = qcid,
    met_scaling = met_scaling
  )
}