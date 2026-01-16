#' exports missing value report as excel file
#'
#' @keywords internal
#' @noRd
export_corr_xlsx <- function(all_corr, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  
  names(all_corr$raw)[names(all_corr$raw) == "col1"] <- "Metabolite 1"
  names(all_corr$raw)[names(all_corr$raw) == "col2"] <- "Metabolite 2"
  names(all_corr$raw)[names(all_corr$raw) == "cor"]  <- "Pearson's r"
  # make column names bold and descriptions with orange background
  bold  <- openxlsx::createStyle(textDecoration = "Bold")
  note  <- openxlsx::createStyle(wrapText = TRUE,
                                 valign = "top",
                                 fgFill = "#f8cbad")
  
  .add_sheet <- function(name) {
    nm <- gsub("[\\[\\]\\*\\?:/\\\\]", "_", name)
    nm <- substr(nm, 1L, 31L)
    openxlsx::addWorksheet(wb, nm)
    nm
  }
  
  shiny::withProgress(message = "Creating metabolite correlation summary...", value = 0, {
    if (all_corr$transformed_included)
      N <- 3
    else
      N <- 2
    
    # sheet 1  Raw Data Metabolite Correlations
    s1 <- .add_sheet("Raw Data")
    txt1 <- paste(
      "Tab Raw Data Correlations. Pearson's r values for all metabolite pairs.",
      "If r = -1 or near -1, the pair have a strong negative linear correlation.",
      "If r = 0, there is no correlation and r values near 0 have weak correlations.",
      "If r = 1 or is close to 1, the pair have a strong positive linear correlation.",
      "n-complete is the number of samples with both metabolite values non-missing."
    )
    openxlsx::writeData(wb,
                        s1,
                        x = txt1,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:10, rows = 1)
    openxlsx::addStyle(
      wb,
      s1,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s1, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s1,
      x = all_corr$raw,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / N, detail = "Saved: Raw Data Metabolite Correlations.")
    
    # Sheet 2: Corrected Data Metabolite Correlations
    names(all_corr$corrected)[names(all_corr$corrected) == "col1"] <- "Metabolite 1"
    names(all_corr$corrected)[names(all_corr$corrected) == "col2"] <- "Metabolite 2"
    names(all_corr$corrected)[names(all_corr$corrected) == "cor"]  <- "Pearson's r"
    s2 <- .add_sheet("Corrected Data")
    txt2 <- paste(
      "Tab Corrected Data Correlations. Pearson's r values for all metabolite pairs.",
      "If r = -1 or near -1, the pair have a strong negative linear correlation.",
      "If r = 0, there is no correlation and r values near 0 have weak correlations.",
      "If r = 1 or is close to 1, the pair have a strong positive linear correlation.",
      "n-complete is the number of samples with both metabolite values non-missing."
    )
    openxlsx::writeData(wb,
                        s2,
                        x = txt2,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s2, cols = 1:10, rows = 1)
    openxlsx::addStyle(
      wb,
      s2,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s2, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s2,
      x = all_corr$corrected,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / N, detail = "Saved: Corrected Date Correlations")
    
    if (all_corr$transformed_included) {
      names(all_corr$transformed)[names(all_corr$transformed) == "col1"] <- "Metabolite 1"
      names(all_corr$transformed)[names(all_corr$transformed) == "col2"] <- "Metabolite 2"
      names(all_corr$transformed)[names(all_corr$transformed) == "cor"]  <- "Pearson's r"
      s3 <- .add_sheet("Transformed Data")
      txt3 <- paste(
        "Tab Transformed and Corrected Data Correlations. Pearson's r values for all metabolite pairs.",
        "If r = -1 or near -1, the pair have a strong negative linear correlation.",
        "If r = 0, there is no correlation and r values near 0 have weak correlations.",
        "If r = 1 or is close to 1, the pair have a strong positive linear correlation.",
        "n-complete is the number of samples with both metabolite values non-missing."
      )
      openxlsx::writeData(wb,
                          s3,
                          x = txt3,
                          startCol = 1,
                          startRow = 1)
      openxlsx::mergeCells(wb, s3, cols = 1:10, rows = 1)
      openxlsx::addStyle(
        wb,
        s3,
        style = note,
        rows = 1,
        cols = 1,
        gridExpand = TRUE
      )
      openxlsx::setRowHeights(wb, s3, rows = 1, heights = 60)
      openxlsx::writeData(
        wb,
        s3,
        x = all_corr$transformed,
        startRow = 3,
        startCol = 1,
        headerStyle = bold
      )
      shiny::incProgress(1 / N, detail = "Saved: Transformed and Corrected Date Correlations")
      
    }
  })
  return(wb)
}