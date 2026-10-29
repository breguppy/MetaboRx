#' Require a package or fail with a clear message
#'
#' @keywords internal
#' @noRd
.require_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "%s is required to %s. Install it with install.packages('%s').",
      pkg, why, pkg
    ), call. = FALSE)
  }
}

#' Metadata columns used throughout the correction workflow
#'
#' @keywords internal
#' @noRd
.correction_metadata_cols <- function() {
  c("sample", "batch", "class", "order")
}

#' Metabolite columns for correction workflow data frames
#'
#' @keywords internal
#' @noRd
.correction_metab_cols <- function(df, metadata_cols = .correction_metadata_cols()) {
  setdiff(names(df), metadata_cols)
}

#' Whether QC RSD filtering is enabled
#'
#' @keywords internal
#' @noRd
.qc_rsd_filter_enabled <- function(p, rsd_cutoff = NULL) {
  if (!is.null(p$remove_qc_rsd_filter)) {
    return(isTRUE(p$remove_qc_rsd_filter))
  }

  if (!is.null(p$post_cor_filter)) {
    return(!isTRUE(p$post_cor_filter))
  }

  length(rsd_cutoff) == 1L && !is.na(rsd_cutoff) && is.finite(rsd_cutoff)
}
#' Final metabolite cleanup used by correction methods
#'
#' @keywords internal
#' @noRd
.cleanup_corrected_metabolites <- function(df, metab_cols) {
  if (length(metab_cols) == 0L) {
    return(df)
  }

  df[metab_cols] <- lapply(df[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    x[is.na(x)] <- if (is.finite(mp)) mp else 0
    x
  })

  df
}

#' Align pre-imputation reference values by unique sample identifier
#'
#' NULL means that the supplied correction input is itself the reference.
#' @keywords internal
#' @noRd
.align_qc_reference <- function(df, reference_df, metab_cols) {
  if (is.null(reference_df)) return(df)
  required <- c("sample", "class", metab_cols)
  if (!is.data.frame(reference_df) ||
      !all(required %in% names(reference_df))) {
    stop("reference_df must contain sample, class, and metabolite columns.")
  }
  if (!all(c("sample", "class") %in% names(df))) {
    stop("Correction input must contain sample and class columns.")
  }
  if (anyNA(df$sample) || anyDuplicated(df$sample) ||
      anyNA(reference_df$sample) || anyDuplicated(reference_df$sample)) {
    stop("Sample identifiers must be nonmissing and unique in both data frames.")
  }
  idx <- match(df$sample, reference_df$sample)
  if (anyNA(idx)) stop("reference_df is missing correction samples.")
  aligned <- reference_df[idx, , drop = FALSE]
  if (!identical(as.character(aligned$class), as.character(df$class))) {
    stop("Reference and correction sample classes must agree.")
  }
  if (!all(vapply(aligned[metab_cols], is.numeric, logical(1)))) {
    stop("Reference metabolite columns must be numeric.")
  }
  aligned
}

#' Normalize, impute, clean, and scale corrected metabolites
#'
#' @keywords internal
#' @noRd
.finalize_qc_correction <- function(
    out,
    raw_df,
    metab_cols,
    qcid,
    met_scaling,
    reference_df = NULL,
    impute_missing = TRUE
) {
  met_scaling <- match.arg(met_scaling, c("one", "original"))
  
  reference_df <- .align_qc_reference(raw_df, reference_df, metab_cols)

  # Retain observed positive QC rows before correction-stage imputation.
  # The raw target uses ALL observed positive QCs, even if a fit fails.
  qc_reference <- vector("list", length(metab_cols))
  names(qc_reference) <- metab_cols
  
  for (metab in metab_cols) {
    raw_y <- as.numeric(raw_df[[metab]])
    corr <- as.numeric(out[[metab]])
    
    corr[!is.finite(corr) | corr < 0] <- NA_real_
    
    observed_y <- as.numeric(reference_df[[metab]])
    observed_qc <- is.finite(observed_y[qcid]) & observed_y[qcid] > 0
    valid_qc <- observed_qc &
      is.finite(corr[qcid]) &
      corr[qcid] > 0
    
    reference_rows <- qcid[valid_qc]
    
    qc_reference[[metab]] <- list(
      rows = reference_rows,
      raw_median = if (any(observed_qc)) {
        stats::median(observed_y[qcid][observed_qc])
      } else {
        NA_real_
      }
    )
    
    # Common scale for imputation, regardless of met_scaling.
    if (length(reference_rows) > 0L) {
      sf <- stats::median(corr[reference_rows])
      corr <- corr / sf
    }
    
    zero_mask <- is.finite(raw_y) & raw_y == 0
    corr[zero_mask] <- 0
    
    out[[metab]] <- corr
  }
  
  # Handle almost entirely missing metabolites before kNN.
  if (isTRUE(impute_missing) && anyNA(out[metab_cols])) {
    for (metab in metab_cols) {
      x <- out[[metab]]
      
      if (mean(is.na(x)) >= 0.99) {
        positive <- x[is.finite(x) & x > 0]
        replacement <- if (length(positive) > 0L) {
          min(positive)
        } else {
          0
        }
        
        x[is.na(x)] <- replacement
        out[[metab]] <- x
      }
    }
    
    if (anyNA(out[metab_cols]) && length(metab_cols) >= 2L) {
      kn <- impute::impute.knn(
        t(as.matrix(out[metab_cols])),
        rowmax = 1,
        colmax = 1,
        maxp = 15000
      )
      
      out[metab_cols] <- as.data.frame(t(kn$data))
    }
  }
  
  out <- .cleanup_corrected_metabolites(out, metab_cols)
  
  # Restore original zeros and apply final scaling after cleanup.
  for (metab in metab_cols) {
    raw_y <- as.numeric(raw_df[[metab]])
    zero_mask <- is.finite(raw_y) & raw_y == 0
    
    x <- as.numeric(out[[metab]])
    x[zero_mask] <- 0
    
    reference <- qc_reference[[metab]]
    reference_rows <- reference$rows
    
    if (length(reference_rows) == 0L) {
      # Entirely zero metabolites have no positive QC scaling target.
      if (any(is.finite(x) & x > 0)) {
        warning(
          sprintf(
            "Could not scale '%s': no valid positive observed QC pairs.",
            metab
          )
        )
      }
      
      out[[metab]] <- x
      next
    }
    
    reference_values <- x[reference_rows]
    
    if (any(!is.finite(reference_values) | reference_values <= 0)) {
      warning(
        sprintf(
          "Could not scale '%s': observed QC reference values became invalid.",
          metab
        )
      )
      
      out[[metab]] <- x
      next
    }
    
    corrected_median <- stats::median(reference_values)
    
    target <- if (met_scaling == "original") {
      reference$raw_median
    } else {
      1
    }
    
    scale_factor <- target / corrected_median
    
    if (!is.finite(scale_factor) || scale_factor <= 0) {
      warning(sprintf("Could not calculate scaling for '%s'.", metab))
      out[[metab]] <- x
      next
    }
    
    scaled <- x * scale_factor
    
    # Avoid introducing overflow through final scaling.
    if (any(is.finite(x) & !is.finite(scaled))) {
      warning(
        sprintf(
          "Could not scale '%s': final scaling would cause numeric overflow.",
          metab
        )
      )
      
      out[[metab]] <- x
      next
    }
    
    scaled[zero_mask] <- 0
    out[[metab]] <- scaled
  }
  
  out
}

# Silence R CMD check notes for NSE/dplyr/ggplot2 variables
utils::globalVariables(c(
  ".", ".data", "Batch", "BatchLabel", "CDF", "Group", "Mean", "Metabolite",
  "PC1", "PC2", "RSD", "RSD_NonQC", "RSD_NonQC_after", "RSD_NonQC_before",
  "RSD_QC", "RSD_QC_after", "RSD_QC_before", "RSD_after", "RSD_before", "SD",
  "Type", "Value", "across", "addStyle", "addWorksheet", "aes", "all_of",
  "batch", "bind_cols", "card_title", "cell_line", "change", "column",
  "conditionalPanel", "createStyle", "createWorkbook", "data_type", "decrease",
  "desc", "distinct", "dose", "downloadButton", "downloadHandler",
  "element_text", "eventReactive", "fill", "fluidRow", "geom_line",
  "geom_point", "ggplot", "ggsave", "guides", "isolate", "labs",
  "left_join", "mean_RSD", "median", "mergeCells", "modifyList", "n", "p",
  "percent", "percent_format", "pivot_longer", "plotOutput", "prcomp", "predict",
  "pull", "radioButtons", "read_excel", "renderPlot", "rsd_QC",
  "rsd_after", "rsd_before", "rsd_nonqc_after", "rsd_nonqc_before",
  "rsd_qc_after", "rsd_qc_before", "saveWorkbook", "scale_color_brewer",
  "scale_color_manual", "scale_x_continuous", "scale_y_continuous", "sd",
  "setColWidths", "setNames", "setRowHeights", "setTxtProgressBar",
  "showNotification", "slice_head", "stopApp", "tagList", "theme",
  "theme_minimal", "tibble", "type", "ungroup", "withProgress", "withSpinner",
  "writeData", "xlim", "xmax", "xmin", "y", "ylim", "zip", "RSD_a", "RSD_b",
  "delta", "flagged", "group_id", "z", "panel", "variable_wrapped", "loading",
  "decision", " ", ":=", "QC_after", "QC_before", "Samples_after",
  "Samples_before", "after", "analysis", "before", "dataset", "group",
  "unit_id"
))
