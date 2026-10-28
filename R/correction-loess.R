#' @keywords internal
#' @noRd
.loess_numerical_warning <- function(messages) {
  if (length(messages) == 0L) {
    return(FALSE)
  }

  any(grepl(
    "pseudoinverse|near singular|reciprocal condition|neighborhood radius|zero-width",
    messages,
    ignore.case = TRUE
  ))
}

#' @keywords internal
#' @noRd
.with_loess_warnings <- function(expr) {
  warning_messages <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(NA_real_, error_message = conditionMessage(e))
    }
  )

  list(value = value, warnings = warning_messages)
}

#' @keywords internal
#' @noRd
.loess_prediction_ok <- function(pred, qc_y, max_fold = 10, small_n_envelope_fold = 1.25) {
  if (!is.numeric(pred) || length(pred) == 0L || any(!is.finite(pred)) || any(pred <= 0)) {
    return(FALSE)
  }

  qc_y <- qc_y[is.finite(qc_y) & qc_y > 0]
  qc_ref <- stats::median(qc_y, na.rm = TRUE)
  if (!is.finite(qc_ref) || qc_ref <= 0) {
    return(FALSE)
  }

  lower <- qc_ref / max_fold
  upper <- qc_ref * max_fold
  if (!all(pred >= lower & pred <= upper)) {
    return(FALSE)
  }

  if (length(qc_y) <= 6L) {
    small_n_lower <- min(qc_y, na.rm = TRUE) / small_n_envelope_fold
    small_n_upper <- max(qc_y, na.rm = TRUE) * small_n_envelope_fold
    if (!all(pred >= small_n_lower & pred <= small_n_upper)) {
      return(FALSE)
    }
  }

  TRUE
}

#' @keywords internal
#' @noRd
.loess_candidate_predict_x <- function(qx, qy, newx, span, degree, family = "symmetric") {
  fit_result <- .with_loess_warnings({
    stats::loess(
      qy ~ qx,
      span = span,
      degree = degree,
      family = family,
      control = stats::loess.control(surface = "direct")
    )
  })

  fit <- fit_result$value
  if (!inherits(fit, "loess")) {
    error_message <- attr(fit, "error_message", exact = TRUE)
    return(list(
      pred = NA_real_,
      ok = FALSE,
      reason = if (is.null(error_message)) "loess_error" else error_message
    ))
  }

  pred_result <- .with_loess_warnings({
    stats::predict(fit, newdata = data.frame(qx = newx))
  })
  qc_result <- .with_loess_warnings({
    stats::predict(fit, newdata = data.frame(qx = qx))
  })

  warning_messages <- c(fit_result$warnings, pred_result$warnings, qc_result$warnings)
  pred <- pred_result$value
  qc_pred <- qc_result$value

  if (.loess_numerical_warning(warning_messages)) {
    return(list(pred = pred, ok = FALSE, reason = "loess_numerical_warning"))
  }
  if (!.loess_prediction_ok(pred, qy) || !.loess_prediction_ok(qc_pred, qy)) {
    return(list(pred = pred, ok = FALSE, reason = "loess_implausible_prediction"))
  }

  list(pred = pred, ok = TRUE, reason = NA_character_, family = family, scale = "raw")
}

#' @keywords internal
#' @noRd
.annotate_loess_prediction <- function(pred, fit_method, fallback_reason = NA_character_, fit_family = NA_character_, fit_scale = NA_character_) {
  attr(pred, "fit_method") <- fit_method
  attr(pred, "fallback_reason") <- fallback_reason
  attr(pred, "fit_family") <- fit_family
  attr(pred, "fit_scale") <- fit_scale
  pred
}

#' @keywords internal
#' @noRd
.safe_loess_predict_x <- function(qc_x, qc_y, newx, span, degree) {
  ok <- is.finite(qc_x) & is.finite(qc_y) & qc_y > 0
  qx <- as.numeric(qc_x[ok])
  qy <- as.numeric(qc_y[ok])

  if (length(qx) < 2L) {
    return(.annotate_loess_prediction(rep(1, length(newx)), "constant", "fewer_than_two_qcs"))
  }

  qc_df <- data.frame(qx = qx, qy = qy)
  qc_df <- stats::aggregate(qy ~ qx, data = qc_df, FUN = stats::median)
  ord <- order(qc_df$qx)
  qx <- as.numeric(qc_df$qx[ord])
  qy <- as.numeric(qc_df$qy[ord])

  if (length(qx) < 2L) {
    return(.annotate_loess_prediction(rep(stats::median(qy, na.rm = TRUE), length(newx)), "constant", "fewer_than_two_distinct_qcs"))
  }

  deg_req <- as.integer(degree)
  deg_req <- max(0L, min(2L, deg_req))
  n <- length(qx)
  spn <- max(as.numeric(span)[1], min(1, 8 / n))
  if (!is.finite(spn) || spn <= 0) {
    spn <- 0.75
  }
  spn <- min(1, spn)

  fallback_reason <- NA_character_
  for (deg in seq.int(deg_req, 0L)) {
    if (length(unique(qx)) < (deg + 1L)) {
      fallback_reason <- "too_few_distinct_qcs_for_degree"
      next
    }

    candidate <- .loess_candidate_predict_x(
      qx = qx,
      qy = qy,
      newx = newx,
      span = spn,
      degree = deg,
      family = "symmetric"
    )
    if (isTRUE(candidate$ok)) {
      method <- sprintf("loess_degree_%d", deg)
      return(.annotate_loess_prediction(candidate$pred, method, fallback_reason, candidate$family, candidate$scale))
    }

    fallback_reason <- candidate$reason
  }

  if (exists(".safe_nw_predict_x", mode = "function")) {
    pred <- .safe_nw_predict_x(qc_x = qx, qc_y = qy, newx = newx, span = span)
    return(.annotate_loess_prediction(pred, "local_constant", fallback_reason, NA_character_, NA_character_))
  }

  pred <- stats::approx(qx, qy, xout = newx, rule = 2)$y
  .annotate_loess_prediction(pred, "linear_interpolation", fallback_reason)
}

loess_correction <- function(df, metab_cols, degree, met_scaling = "one", span = 0.75, min_qc = 3) {
  df <- df[order(df$order), , drop = FALSE]
  if (!(identical(df$class[1], "QC") && identical(df$class[nrow(df)], "QC"))) {
    stop("First and last samples must be QCs.")
  }

  qcid <- which(df$class == "QC")
  if (length(qcid) < min_qc) stop(sprintf("Need at least %d QC rows for LOESS.", min_qc))

  x_all <- suppressWarnings(as.numeric(df$order))
  if (any(!is.finite(x_all))) stop("order must be numeric and finite after sorting.")

  out <- df
  diagnostics <- data.frame(
    metabolite = metab_cols,
    fit_method = rep(NA_character_, length(metab_cols)),
    fallback_reason = rep(NA_character_, length(metab_cols)),
    fit_family = rep(NA_character_, length(metab_cols)),
    fit_scale = rep(NA_character_, length(metab_cols)),
    stringsAsFactors = FALSE
  )

  for (metab in metab_cols) {
    zero_mask <- is.finite(df[[metab]]) & df[[metab]] == 0
    qc_y <- df[[metab]][qcid]

    if (all(qc_y <= 0, na.rm = TRUE)) {
      out[[metab]] <- 0
      diag_idx <- match(metab, diagnostics$metabolite)
      diagnostics$fit_method[diag_idx] <- "all_nonpositive_qc"
      diagnostics$fallback_reason[diag_idx] <- "all_nonpositive_qc"
      diagnostics$fit_family[diag_idx] <- NA_character_
      diagnostics$fit_scale[diag_idx] <- NA_character_
      next
    }

    pred <- .safe_loess_predict_x(
      qc_x = x_all[qcid],
      qc_y = qc_y,
      newx = x_all,
      span = span,
      degree = degree
    )

    diag_idx <- match(metab, diagnostics$metabolite)
    diagnostics$fit_method[diag_idx] <- attr(pred, "fit_method", exact = TRUE)
    diagnostics$fallback_reason[diag_idx] <- attr(pred, "fallback_reason", exact = TRUE)
    diagnostics$fit_family[diag_idx] <- attr(pred, "fit_family", exact = TRUE)
    diagnostics$fit_scale[diag_idx] <- attr(pred, "fit_scale", exact = TRUE)

    pred[!is.finite(pred) | pred <= 0] <- NA_real_
    corr <- as.numeric(df[[metab]]) / pred

    if(met_scaling == "one") {
      sf <- stats::median(corr[qcid], na.rm = TRUE)
      if (is.finite(sf) && sf > 0) {
        corr <- corr / sf
      }
    } else if(met_scaling == "original"){
      valid_qc <- qc_y[is.finite(qc_y) & qc_y > 0]
      qc_reference <- stats::median(valid_qc)
      
      sf <- stats::median(corr[qcid], na.rm = TRUE)
      if (is.finite(sf) && sf > 0 && is.finite(qc_reference)) {
        corr <- (corr / sf) * qc_reference
      }
    }

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

  out <- .cleanup_corrected_metabolites(out, metab_cols)
  attr(out, "loess_diagnostics") <- diagnostics
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
  batch_ids <- unique(df$batch)
  rows_by_batch <- lapply(batch_ids, function(b) which(df$batch == b))
  diagnostics <- data.frame(
    metabolite = character(0),
    batch = character(0),
    fit_method = character(0),
    fallback_reason = character(0),
    fit_family = character(0),
    fit_scale = character(0),
    stringsAsFactors = FALSE
  )

  for (metab in metab_cols) {
    # preserve original exact zeros (common for missing/below-LOD encoding)
    zero_mask <- is.finite(df[[metab]]) & df[[metab]] == 0

    for (i in seq_along(batch_ids)) {
      b <- batch_ids[[i]]
      b_idx <- rows_by_batch[[i]]
      sub <- df[b_idx, , drop = FALSE]

      if (!(identical(sub$class[1], "QC") && identical(sub$class[nrow(sub)], "QC"))) {
        stop(sprintf("Batch '%s' must start and end with QC.", b))
      }

      qcid <- which(sub$class == "QC")
      if (length(qcid) < min_qc) {
        warning(sprintf(
          "Skipping batch '%s' for '%s': only %d QC rows (< %d).",
          b, metab, length(qcid), min_qc
        ))
        diagnostics <- rbind(
          diagnostics,
          data.frame(
            metabolite = metab,
            batch = as.character(b),
            fit_method = "skipped_too_few_qc",
            fallback_reason = "too_few_qc",
            fit_family = NA_character_,
            fit_scale = NA_character_,
            stringsAsFactors = FALSE
          )
        )
        next
      }

      qc_y <- sub[[metab]][qcid]
      if (all(qc_y <= 0, na.rm = TRUE)) {
        out[[metab]][b_idx] <- 0
        diagnostics <- rbind(
          diagnostics,
          data.frame(
            metabolite = metab,
            batch = as.character(b),
            fit_method = "all_nonpositive_qc",
            fallback_reason = "all_nonpositive_qc",
            fit_family = NA_character_,
            fit_scale = NA_character_,
            stringsAsFactors = FALSE
          )
        )
        next
      }

      x_sub <- suppressWarnings(as.numeric(sub$order))
      if (any(!is.finite(x_sub))) {
        stop(sprintf("Non-numeric or non-finite `order` detected in batch '%s'.", b))
      }

      pred <- .safe_loess_predict_x(
        qc_x = x_sub[qcid],
        qc_y = qc_y,
        newx = x_sub,
        span = span,
        degree = degree
      )

      diagnostics <- rbind(
        diagnostics,
        data.frame(
          metabolite = metab,
          batch = as.character(b),
          fit_method = attr(pred, "fit_method", exact = TRUE),
          fallback_reason = attr(pred, "fallback_reason", exact = TRUE),
          fit_family = attr(pred, "fit_family", exact = TRUE),
          fit_scale = attr(pred, "fit_scale", exact = TRUE),
          stringsAsFactors = FALSE
        )
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
        maxp = 15000
      )
      out[metab_cols] <- as.data.frame(t(kn$data))
    }
  }

  # Final cleanup: enforce non-negative finite values, fill remaining NA with
  # smallest positive, else 0.
  out <- .cleanup_corrected_metabolites(out, metab_cols)
  attr(out, "loess_diagnostics") <- diagnostics
  out
}
