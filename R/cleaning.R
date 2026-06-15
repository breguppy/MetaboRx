#' Clean data and prep for correction
#'
#' This function:
#' 1) Standardizes required metadata column names
#' 2) Removes additional metadata columns from the cleaned data frame and saves
#'    them in `meta_df`
#' 3) Removes non-numerical metabolite columns and columns containing all missing
#'    values or zeros for QC samples
#' 4) Removes HP rows and separates blank/PB rows
#' 5) Counts non-numerical and exact zeros that are replaced as missing values
#' 6) Makes sure injection order starts and ends with a QC sample
#' 7) Finds duplicate metabolite columns
#'
#' @param df A data frame containing metabolite columns of the raw data.
#' @param sample User-selected column name containing sample names.
#' @param batch User-selected column name containing batch information.
#'   If a batch column is not provided one will be made.
#' @param class User-selected column name containing sample classification
#'   groups.
#' @param order User-selected column name containing sample injection order.
#' @param withheld_cols A vector of user-selected additional metadata columns.
#'
#' @return A list containing cleaned data, metadata, removed columns, replacement
#'   counts, separated blanks, and duplicate metabolite information.
#'
#' @keywords internal
#' @noRd
clean_data <- function(df,
                       sample,
                       batch,
                       class,
                       order,
                       withheld_cols) {
  col_names <- names(df)

  duplicate_col_names <- unique(col_names[duplicated(col_names)])

  if (length(duplicate_col_names) > 0L) {
    df <- repair_duplicate_column_names(df)
  }

  # 1. Standardize required metadata column names ------------------------------
  if (!(batch %in% names(df))) {
    df[["batch"]] <- "batch1"
  }

  df_names <- names(df)
  df_names[df_names == sample] <- "sample"
  df_names[df_names == batch] <- "batch"
  df_names[df_names == class] <- "class"
  df_names[df_names == order] <- "order"
  names(df) <- df_names

  df[["class"]] <- as.character(df[["class"]])
  df[["order"]] <- as.numeric(df[["order"]])

  df <- df[order(df[["order"]]), , drop = FALSE]

  class_chr <- trimws(as.character(df[["class"]]))
  class_chr[is.na(class_chr) | class_chr == ""] <- "QC"

  qc_idx_norm <- tolower(class_chr) == "qc"
  class_chr[qc_idx_norm] <- "QC"

  df[["class"]] <- class_chr

  # 2. Remove additional metadata columns and save them ------------------------
  meta_columns <- intersect(
    c("sample", "batch", "class", "order", withheld_cols),
    names(df)
  )

  meta_df <- df[, meta_columns, drop = FALSE]

  df <- df[, setdiff(names(df), withheld_cols), drop = FALSE]

  # 3. Remove invalid metabolite columns and save names ------------------------
  metab_candidates <- setdiff(names(df), c("sample", "batch", "class", "order"))

  non_numeric_cols <- metab_candidates[
    vapply(
      df[, metab_candidates, drop = FALSE],
      FUN = function(col) {
        vals <- col[!is.na(col)]
        all(is.na(suppressWarnings(as.numeric(vals))))
      },
      FUN.VALUE = logical(1L)
    )
  ]

  metab_numeric <- setdiff(metab_candidates, non_numeric_cols)

  is_all_missing_or_zero <- function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    all(is.na(x_num) | x_num == 0)
  }

  qc_idx_for_filter <- trimws(as.character(df[["class"]])) == "QC"

  all_missing_zero_qc_cols <- character(0)

  if (length(metab_numeric) > 0L && any(qc_idx_for_filter)) {
    all_missing_zero_qc_cols <- metab_numeric[
      vapply(
        df[qc_idx_for_filter, metab_numeric, drop = FALSE],
        FUN = is_all_missing_or_zero,
        FUN.VALUE = logical(1L)
      )
    ]
  }

  df <- df[
    ,
    !(names(df) %in% c(non_numeric_cols, all_missing_zero_qc_cols)),
    drop = FALSE
  ]

  # 4. Remove HP rows and separate blanks/PBs ----------------------------------
  is_hp <- toupper(trimws(df[["class"]])) == "HP"

  if (any(is_hp, na.rm = TRUE)) {
    df <- df[!is_hp, , drop = FALSE]
  }

  df <- df[order(df[["order"]]), , drop = FALSE]
  
  blank_like_labels <- c("blank", "pb", "processing blank")
  class_chr <- trimws(as.character(df[["class"]]))
  is_blank_like <- tolower(class_chr) %in% blank_like_labels

  blank_df <- df[0, , drop = FALSE]

  if (any(is_blank_like, na.rm = TRUE)) {
    blank_df <- df[is_blank_like, , drop = FALSE]
    df <- df[!is_blank_like, , drop = FALSE]
    df <- df[order(df[["order"]]), , drop = FALSE]
  }

  # Recalculate metabolite columns after row filtering.
  metab <- setdiff(names(df), c("sample", "batch", "class", "order"))

  # 5. Count non-numerical and zeros replaced with NA --------------------------
  if (length(metab) > 0L) {
    original_metab <- df[, metab, drop = FALSE]
    numeric_metab <- lapply(original_metab, function(col) {
      suppressWarnings(as.numeric(col))
    })

    non_numeric_replaced <- vapply(
      seq_along(metab),
      function(i) {
        sum(is.na(numeric_metab[[i]]) & !is.na(original_metab[[i]]))
      },
      integer(1L)
    )

    zero_replaced <- vapply(
      numeric_metab,
      function(col) {
        sum(col == 0, na.rm = TRUE)
      },
      integer(1L)
    )

    numeric_metab <- lapply(numeric_metab, function(col) {
      zero_idx <- !is.na(col) & col == 0
      col[zero_idx] <- NA_real_
      col
    })

    df[, metab] <- numeric_metab
  } else {
    non_numeric_replaced <- integer(0)
    zero_replaced <- integer(0)
  }

  repl <- tibble::tibble(
    metabolite = metab,
    non_numeric_replaced = non_numeric_replaced,
    zero_replaced = zero_replaced
  )

  # 6. Make sure data starts and ends with a QC --------------------------------
  if (nrow(df) == 0L) {
    stop("No non-blank/non-HP rows remain after preprocessing.")
  }

  if (df[["class"]][1] != "QC") {
    stop("Data sorted by injection order must begin with a QC sample.")
  }

  if (df[["class"]][nrow(df)] != "QC") {
    stop("Data sorted by injection order must end with a QC sample.")
  }

  # 7. Find equal columns ------------------------------------------------------
  duplicate_mets <- find_equal_metabolite_cols(df, metab, tolerance = 1e-3)

  list(
    df = df,
    meta_df = meta_df,
    replacement_counts = repl,
    withheld_cols = withheld_cols,
    non_numeric_cols = non_numeric_cols,
    all_missing_zero_qc_cols = all_missing_zero_qc_cols,
    duplicate_mets = duplicate_mets,
    duplicate_col_names = duplicate_col_names,
    blank_df = blank_df
  )
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

  col_values <- lapply(cols, function(col) df[[col]])
  names(col_values) <- cols
  col_present <- lapply(col_values, function(col) !is.na(col))
  results <- list()
  result_idx <- 0L

  for (i in seq_len(length(cols) - 1L)) {
    x <- col_values[[i]]
    x_present <- col_present[[i]]

    for (j in seq.int(i + 1L, length(cols))) {
      idx <- x_present & col_present[[j]]

      if (!any(idx)) {
        next
      }

      if (isTRUE(all.equal(x[idx], col_values[[j]][idx], ...))) {
        result_idx <- result_idx + 1L
        results[[result_idx]] <- data.frame(
          col1 = cols[[i]],
          col2 = cols[[j]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(results)) {
    return(data.frame(col1 = character(0), col2 = character(0)))
  }

  do.call(rbind, results)
}

#' Repair duplicated column names safely
#'
#' Keeps the first occurrence of each duplicated column name unchanged and
#' appends "_1", "_2", etc. to subsequent duplicates. Generated names are
#' guaranteed not to collide with existing column names.
#'
#' @param df A data frame.
#'
#' @return A data frame with unique column names.
#'
#' @keywords internal
#' @noRd
repair_duplicate_column_names <- function(df) {
  col_names <- names(df)
  
  if (length(col_names) == 0L) {
    return(df)
  }
  
  repaired_names <- character(length(col_names))
  used_names <- character(0)
  duplicate_counts <- integer(0)
  names(duplicate_counts) <- character(0)
  
  for (i in seq_along(col_names)) {
    current_name <- col_names[[i]]
    
    if (!current_name %in% names(duplicate_counts)) {
      duplicate_counts[[current_name]] <- 0L
    }
    
    if (!current_name %in% used_names) {
      repaired_names[[i]] <- current_name
      used_names <- c(used_names, current_name)
      next
    }
    
    repeat {
      duplicate_counts[[current_name]] <- duplicate_counts[[current_name]] + 1L
      candidate_name <- paste0(current_name, "_", duplicate_counts[[current_name]])
      
      if (!candidate_name %in% used_names && !candidate_name %in% col_names) {
        repaired_names[[i]] <- candidate_name
        used_names <- c(used_names, candidate_name)
        break
      }
    }
  }
  
  names(df) <- repaired_names
  
  df
}
