#' Compute pairwise correlations for metabolite columns, ignoring NAs/infinite values
#'
#' @param df           Data frame containing metabolite columns (and possibly metadata).
#' @param cols         Optional character vector of column names to consider. If NULL, uses all columns in `df`.
#' @param method       Correlation method passed to stats::cor() (e.g., "pearson", "spearman", "kendall").
#' @param min_complete Minimum number of paired finite observations required to compute a correlation.
#'
#' @return A data.frame with columns:
#'   - col1, col2       : column name pairs
#'   - cor              : correlation value
#'   - n_complete       : number of paired finite observations used
#' sorted by descending absolute correlation.
#'
#' @keywords internal
#' @noRd
compute_pairwise_metabolite_correlations <- function(df,
                                                     cols = NULL,
                                                     method = "pearson",
                                                     min_complete = 3L) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(method), length(method) == 1L)
  stopifnot(is.numeric(min_complete), length(min_complete) == 1L, min_complete >= 1L)

  if (is.null(cols)) cols <- names(df)
  cols <- intersect(cols, names(df))

  # Restrict to numeric columns
  is_num <- vapply(df[cols], is.numeric, logical(1L))
  cols <- cols[is_num]

  empty <- data.frame(
    col1 = character(0),
    col2 = character(0),
    cor = numeric(0),
    n_complete = integer(0),
    stringsAsFactors = FALSE
  )

  if (length(cols) < 2L) {
    return(empty)
  }

  pair_idx <- utils::combn(seq_along(cols), 2L)
  data_mat <- as.matrix(df[, cols, drop = FALSE])

  n_pairs <- ncol(pair_idx)
  col1 <- character(n_pairs)
  col2 <- character(n_pairs)
  cor_values <- numeric(n_pairs)
  n_complete_values <- integer(n_pairs)
  keep <- logical(n_pairs)

  for (i in seq_len(n_pairs)) {
    i1 <- pair_idx[1L, i]
    i2 <- pair_idx[2L, i]
    c1 <- cols[[i1]]
    c2 <- cols[[i2]]

    x <- data_mat[, i1]
    y <- data_mat[, i2]

    idx <- is.finite(x) & is.finite(y)
    n_complete <- sum(idx)

    if (n_complete < min_complete) next

    cor_val <- suppressWarnings(stats::cor(x[idx], y[idx], method = method))
    if (is.na(cor_val)) next

    col1[[i]] <- c1
    col2[[i]] <- c2
    cor_values[[i]] <- unname(cor_val)
    n_complete_values[[i]] <- as.integer(n_complete)
    keep[i] <- TRUE
  }

  if (!any(keep)) {
    return(empty)
  }

  out <- data.frame(
    col1 = col1[keep],
    col2 = col2[keep],
    cor = cor_values[keep],
    n_complete = n_complete_values[keep],
    stringsAsFactors = FALSE
  )
  out[order(-abs(out$cor)), , drop = FALSE]
}

#' Filter correlation pairs by correlation range
#'
#' @param corr_df A data.frame produced by compute_pairwise_metabolite_correlations().
#' @param range   Numeric vector of length 2: c(min, max). Inclusive.
#'
#' @return A data.frame subset of `corr_df` where cor is in [min, max],
#' sorted by descending absolute correlation.
#'
#' @keywords internal
#' @noRd
filter_correlation_pairs_by_range <- function(corr_df, range) {
  stopifnot(is.data.frame(corr_df))
  stopifnot(is.numeric(range), length(range) == 2L, all(is.finite(range)))

  if (!all(c("col1", "col2", "cor", "n_complete") %in% names(corr_df))) {
    stop("corr_df must have columns: col1, col2, cor, n_complete.")
  }

  rmin <- min(range)
  rmax <- max(range)

  out <- corr_df[corr_df$cor >= rmin & corr_df$cor <= rmax, , drop = FALSE]
  out[order(-abs(out$cor)), , drop = FALSE]
}
