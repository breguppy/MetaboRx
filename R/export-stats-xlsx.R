# ------------------------------------------------------------------------------
# Excel export helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.write_rsd_sheet <- function(wb,
                             sheet_name,
                             x,
                             note_text,
                             bold_style,
                             note_style) {
  openxlsx::addWorksheet(wb, sheet_name)
  
  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = note_text,
    startCol = 1,
    startRow = 1
  )
  
  openxlsx::mergeCells(wb, sheet_name, cols = 1:12, rows = 1)
  
  openxlsx::addStyle(
    wb,
    sheet_name,
    style = note_style,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  
  openxlsx::setRowHeights(wb, sheet_name, rows = 1, heights = 40)
  
  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = x,
    startRow = 3,
    headerStyle = bold_style
  )
}

#' exports corrected data as excel file
#'
#' @keywords internal
#' @noRd
export_stats_xlsx <- function(p, d, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  
  wb <- openxlsx::createWorkbook()
  
  bold <- openxlsx::createStyle(textDecoration = "Bold")
  note <- openxlsx::createStyle(
    wrapText = TRUE,
    valign = "top",
    fgFill = "#f8cbad"
  )
  
  rsd_data <- .get_rsd_data_after(p$rsd_compare, p, d)
  rsd_results <- .build_rsd_results(rsd_data$df_before, rsd_data$df_after)
  
  txt_met <- paste(
    "RSD values are computed for each metabolite.",
    "For each metabolite, RSD = 100% * (standard deviation) / (mean)."
  )
  
  txt_class <- paste(
    "RSD values are computed for each metabolite and class.",
    "For each metabolite, RSD = 100% * (class standard deviation) / (class mean)."
  )
  
  df_compare_met <- rsd_results$metabolite$compare |>
    dplyr::select(Metabolite, Type, before, after, delta) |>
    dplyr::rename(
      RSD_before = before,
      RSD_after = after,
      delta_RSD = delta
    )
  
  df_compare_class <- rsd_results$class_metabolite$compare |>
    dplyr::select(Metabolite, class, Type, before, after, delta) |>
    dplyr::rename(
      RSD_before = before,
      RSD_after = after,
      delta_RSD = delta
    )
  
  shiny::withProgress(message = "Creating rsd_stats_*today's_date*.xlsx...", value = 0, {
    .write_rsd_sheet(
      wb = wb,
      sheet_name = "Raw Metabolite RSD",
      x = rsd_results$metabolite$before,
      note_text = txt_met,
      bold_style = bold,
      note_style = note
    )
    shiny::incProgress(1 / 6, detail = "Saved: Raw Metabolite RSD")
    
    .write_rsd_sheet(
      wb = wb,
      sheet_name = paste(rsd_data$sheet_label, "Metabolite RSD"),
      x = rsd_results$metabolite$after,
      note_text = txt_met,
      bold_style = bold,
      note_style = note
    )
    shiny::incProgress(1 / 6, detail = "Saved: Corrected Metabolite RSD")
    
    .write_rsd_sheet(
      wb = wb,
      sheet_name = "Metabolite RSD Comparison",
      x = df_compare_met,
      note_text = paste(
        txt_met,
        "The change in RSD (delta = after - before) values will be negative for a decrease in RSD, positive for an increase in RSD, and zero for no change in RSD."
      ),
      bold_style = bold,
      note_style = note
    )
    shiny::incProgress(1 / 6, detail = "Saved: Metabolite RSD Comparison")
    
    .write_rsd_sheet(
      wb = wb,
      sheet_name = "Raw Class RSD",
      x = rsd_results$class_metabolite$before,
      note_text = txt_class,
      bold_style = bold,
      note_style = note
    )
    shiny::incProgress(1 / 6, detail = "Saved: Raw Class RSD")
    
    .write_rsd_sheet(
      wb = wb,
      sheet_name = paste(rsd_data$sheet_label, "Class RSD"),
      x = rsd_results$class_metabolite$after,
      note_text = txt_class,
      bold_style = bold,
      note_style = note
    )
    shiny::incProgress(1 / 6, detail = "Saved: Corrected Class RSD")
    
    .write_rsd_sheet(
      wb = wb,
      sheet_name = "Class RSD Comparison",
      x = df_compare_class,
      note_text = paste(
        txt_class,
        "The change in RSD (delta = after - before) values will be negative for a decrease in RSD, positive for an increase in RSD, and zero for no change in RSD."
      ),
      bold_style = bold,
      note_style = note
    )
    shiny::incProgress(1 / 6, detail = "Saved: Class RSD Comparison")
    
    if (!is.null(file)) {
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      return(normalizePath(file, winslash = "/"))
    }
  })
  
  wb
}