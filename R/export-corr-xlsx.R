#' exports missing value report as excel file
#'
#' @keywords internal
#' @noRd
export_corr_xlsx <- function(all_corr, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  styles <- .xlsx_export_styles()

  corr_export <- list(
    raw = .rename_correlation_export_cols(all_corr$raw),
    corrected = .rename_correlation_export_cols(all_corr$corrected),
    transformed = if (isTRUE(all_corr$transformed_included)) {
      .rename_correlation_export_cols(all_corr$transformed)
    } else {
      NULL
    }
  )

  correlation_description <- function(label) {
    paste(
      sprintf("Tab %s Correlations. Pearson's r values for all metabolite pairs.", label),
      "If r = -1 or near -1, the pair have a strong negative linear correlation.",
      "If r = 0, there is no correlation and r values near 0 have weak correlations.",
      "If r = 1 or is close to 1, the pair have a strong positive linear correlation.",
      "n-complete is the number of samples with both metabolite values non-missing."
    )
  }

  write_correlation_sheet <- function(sheet_name, description, x) {
    sheet <- .xlsx_add_sheet(wb, sheet_name)

    openxlsx::writeData(
      wb,
      sheet,
      x = description,
      startCol = 1,
      startRow = 1
    )
    openxlsx::mergeCells(wb, sheet, cols = 1:10, rows = 1)
    openxlsx::addStyle(
      wb,
      sheet,
      style = styles$note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 60)

    openxlsx::writeData(
      wb,
      sheet,
      x = data.frame(t(names(x)), check.names = FALSE),
      startRow = 3,
      startCol = 1,
      colNames = FALSE
    )
    openxlsx::addStyle(
      wb,
      sheet,
      style = styles$bold,
      rows = 3,
      cols = seq_along(x),
      gridExpand = TRUE
    )

    openxlsx::writeData(
      wb,
      sheet,
      x = x,
      startRow = 4,
      startCol = 1,
      colNames = FALSE
    )

    invisible(sheet)
  }

  shiny::withProgress(message = "Creating metabolite correlation summary...", value = 0, {
    if (all_corr$transformed_included) {
      N <- 3
    } else {
      N <- 2
    }

    # sheet 1  Raw Data Metabolite Correlations
    write_correlation_sheet(
      sheet_name = "Raw Data",
      description = correlation_description("Raw Data"),
      x = corr_export$raw
    )
    shiny::incProgress(1 / N, detail = "Saved: Raw Data Metabolite Correlations.")

    # Sheet 2: Corrected Data Metabolite Correlations
    write_correlation_sheet(
      sheet_name = "Corrected Data",
      description = correlation_description("Corrected Data"),
      x = corr_export$corrected
    )
    shiny::incProgress(1 / N, detail = "Saved: Corrected Date Correlations")

    if (all_corr$transformed_included) {
      write_correlation_sheet(
        sheet_name = "Transformed Data",
        description = correlation_description("Transformed and Corrected Data"),
        x = corr_export$transformed
      )
      shiny::incProgress(1 / N, detail = "Saved: Transformed and Corrected Date Correlations")
    }
  })

  .xlsx_save_or_return(wb, file)
}
