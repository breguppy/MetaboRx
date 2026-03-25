#' Clean data and track replacements
#'
#' @keywords internal
#' @noRd
clean_data <- function(df,
                       sample,
                       batch,
                       class,
                       order,
                       withheld_cols) {
  if (!(batch %in% colnames(df))) {
    df$batch <- "batch1"
  }
  
  # remove columns that will not be corrected and are not metadata
  df <- df[, setdiff(names(df), withheld_cols), drop = FALSE]
  
  # Rename metadata columns
  names(df)[names(df) == sample] <- "sample"
  names(df)[names(df) == batch]  <- "batch"
  names(df)[names(df) == class]  <- "class"
  names(df)[names(df) == order]  <- "order"
  
  # make class/order consistent
  df$class <- as.character(df$class)
  df$order <- as.numeric(df$order)
  
  # get names of non-numeric columns that are not metadata columns
  non_numeric_cols <- names(df)[vapply(
    df,
    function(col) {
      vals <- col[!is.na(col)]
      all(is.na(suppressWarnings(as.numeric(vals))))
    },
    logical(1L)
  )]
  non_numeric_cols <- setdiff(non_numeric_cols, c("sample", "batch", "class", "order"))
  
  df <- df[, !(names(df) %in% non_numeric_cols), drop = FALSE]
  
  # make sure data is in injection order
  df <- df[order(df$order), , drop = FALSE]
  metab <- setdiff(names(df), c("sample", "batch", "class", "order"))
  
  # ---- Normalize class labels ----
  class_chr <- trimws(as.character(df$class))
  class_chr[is.na(class_chr) | class_chr == ""] <- "QC"
  
  qc_idx_norm <- tolower(class_chr) %in% c("qc")
  class_chr[qc_idx_norm] <- "QC"
  
  df$class <- class_chr
  
  # ---- Remove HP rows BEFORE building blank_df ----
  is_hp <- toupper(trimws(df$class)) == "HP"
  if (any(is_hp, na.rm = TRUE)) {
    df <- df[!is_hp, , drop = FALSE]
  }
  
  # keep df in injection order after HP removal
  df <- df[order(df$order), , drop = FALSE]
  
  # ---- Define blank-like rows and remove them before replacement counting ----
  blank_like_labels <- c("blank", "pb", "processing blank")
  class_chr <- trimws(as.character(df$class))
  is_blank_like <- tolower(class_chr) %in% blank_like_labels
  
  blank_df <- df[0, , drop = FALSE]
  below_blank_threshold <- character(0)
  below_blank_threshold_ex_ISTD <- character(0)
  
  if (any(is_blank_like, na.rm = TRUE)) {
    blank_df <- df[is_blank_like, , drop = FALSE]
    df <- df[!is_blank_like, , drop = FALSE]
    
    # keep df in injection order after blank/PB removal
    df <- df[order(df$order), , drop = FALSE]
    
    qc_idx <- trimws(as.character(df$class)) == "QC"
    if (!any(qc_idx)) {
      stop("No QC samples remain after removing blanks/PBs; cannot compute blank threshold.")
    }
    
    blank_means <- vapply(
      blank_df[, metab, drop = FALSE],
      function(x) mean(suppressWarnings(as.numeric(x)), na.rm = TRUE),
      numeric(1L)
    )
    
    qc_means <- vapply(
      df[qc_idx, metab, drop = FALSE],
      function(x) mean(suppressWarnings(as.numeric(x)), na.rm = TRUE),
      numeric(1L)
    )
    
    # flag metabolites where QC mean < 3 * blank/PB mean
    eligible <- is.finite(blank_means) & !is.na(blank_means) & (blank_means > 0)
    below_blank_threshold <- names(qc_means)[eligible & (qc_means < (3 * blank_means))]
  }
  
  # ---- Exclude ISTD / ITSD metabolites from threshold flag ----
  below_blank_threshold_ex_ISTD <- below_blank_threshold[
    !grepl("ISTD|ITSD", below_blank_threshold, ignore.case = TRUE)
  ]
  
  # ---- Replacement counting happens AFTER HP + blank/PB removal ----
  repl <- tibble::tibble(
    metabolite = metab,
    non_numeric_replaced = 0L,
    zero_replaced = 0L
  )
  
  for (i in seq_along(metab)) {
    col <- metab[i]
    orig <- df[[col]]
    num <- suppressWarnings(as.numeric(orig))
    
    cnt1 <- sum(is.na(num) & !is.na(orig))
    cnt2 <- sum(num == 0, na.rm = TRUE)
    
    num[num == 0] <- NA
    df[[col]] <- num
    
    repl$non_numeric_replaced[i] <- cnt1
    repl$zero_replaced[i] <- cnt2
  }
  
  # make sure data starts and ends with a QC
  if (nrow(df) == 0L) {
    stop("No non-blank/non-HP rows remain after preprocessing.")
  }
  if (df$class[1] != "QC") {
    stop("Data sorted by injection order must begin with a QC sample.")
  }
  if (df$class[nrow(df)] != "QC") {
    stop("Data sorted by injection order must end with a QC sample.")
  }
  
  # equal columns
  duplicate_mets <- find_equal_metabolite_cols(df, metab, tolerance = 1e-3)
  
  return(list(
    df = df,
    replacement_counts = repl,
    withheld_cols = withheld_cols,
    non_numeric_cols = non_numeric_cols,
    duplicate_mets = duplicate_mets,
    blank_df = blank_df,
    below_blank_threshold = below_blank_threshold,
    below_blank_threshold_ex_ISTD = below_blank_threshold_ex_ISTD
  ))
}

#' Find (nearly) equal columns ignoring NAs
#'
#' @param df   A data frame containing metabolite columns.
#' @param cols Optional character vector of column names to check.
#'             If NULL, all columns in `df` are used.
#' @param ...  Additional arguments passed to `all.equal()`
#'             (e.g. tolerance = 1e-8).
#'
#' @return A data.frame with columns `col1` and `col2` listing ordered pairs
#'         of columns that are equal (per `all.equal`) on all rows where
#'         both columns are non-NA. If no pairs are found, returns an empty
#'         data.frame with the same columns.
#' @keywords internal
#' @noRd
find_equal_metabolite_cols <- function(df, cols = NULL, ...) {
  if (is.null(cols)) {
    cols <- names(df)
  }
  
  cols <- intersect(cols, names(df))
  if (length(cols) < 2L) {
    return(data.frame(col1 = character(0), col2 = character(0)))
  }
  
  pairs <- utils::combn(cols, 2L, simplify = FALSE)
  
  results <- vector("list", length(pairs))
  keep <- logical(length(pairs))
  
  for (i in seq_along(pairs)) {
    c1 <- pairs[[i]][1L]
    c2 <- pairs[[i]][2L]
    
    x <- df[[c1]]
    y <- df[[c2]]
    
    idx <- !is.na(x) & !is.na(y)
    if (!any(idx)) {
      next
    }
    
    if (isTRUE(all.equal(x[idx], y[idx], ...))) {
      results[[i]] <- data.frame(
        col1 = c1,
        col2 = c2,
        stringsAsFactors = FALSE
      )
      keep[i] <- TRUE
    }
  }
  
  if (!any(keep)) {
    return(data.frame(col1 = character(0), col2 = character(0)))
  }
  
  do.call(rbind, results[keep])
}