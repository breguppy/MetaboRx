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
#' Corrects each metabolite by dividing observed values by a
#' Nadaraya-Watson-smoothed QC trend over injection order.
#'
#' @param df data.frame
#'   Data frame containing metadata columns `class` and `order`, plus
#'   metabolite columns.
#' @param metab_cols character vector
#'   Metabolite column names.
#' @param span numeric scalar
#'   Relative smoothing span in (0, 1]. Larger values give smoother fits.
#' @param min_qc integer scalar
#'   Minimum number of QC rows required.
#' @param kernel character scalar
#'   Kernel to use. One of "gaussian" or "tricube".
#'
#' @return data.frame
#'   Data frame with corrected metabolite columns.
#'
#' @keywords internal
#' @noRd
nw_correction <- function(df, metab_cols, span = 0.75, min_qc = 3, kernel = "gaussian") {
  df <- df[order(df$order), , drop = FALSE]

  if (!(identical(df$class[1], "QC") && identical(df$class[nrow(df)], "QC"))) {
    stop("First and last samples must be QCs.")
  }

  qcid <- which(df$class == "QC")
  if (length(qcid) < min_qc) {
    stop(sprintf("Need at least %d QC rows for local constant correction.", min_qc))
  }

  x_all <- suppressWarnings(as.numeric(df$order))
  if (any(!is.finite(x_all))) {
    stop("order must be numeric and finite after sorting.")
  }

  out <- df

  for (metab in metab_cols) {
    zero_mask <- is.finite(df[[metab]]) & df[[metab]] == 0
    qc_y <- df[[metab]][qcid]

    if (all(qc_y <= 0, na.rm = TRUE)) {
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
    corr <- as.numeric(df[[metab]]) / pred

    sf <- stats::median(corr[qcid], na.rm = TRUE)
    if (is.finite(sf) && sf > 0) {
      corr <- corr / sf
    }

    corr[!is.finite(corr) | corr < 0] <- NA_real_
    out[[metab]] <- corr
    out[[metab]][zero_mask] <- 0
  }

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
        maxp = 15000
      )
      out[metab_cols] <- as.data.frame(t(kn$data))
    }
  }

  .cleanup_corrected_metabolites(out, metab_cols)
}
