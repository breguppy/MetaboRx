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
