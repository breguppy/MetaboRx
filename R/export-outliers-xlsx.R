#' Export outliers excel file
#' @keywords internal
#' @noRd
export_outliers_xlsx <- function(p, d, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  # make column names bold and descriptions with orange backgroud
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
  
  df <- d$filtered_corrected$df_no_mv

  res <- detect_hotelling_nonqc_dual_z(df, p)
  
  outlier_samples <- unique(res$data$sample[res$data$is_outlier_sample])
  pc_loadings <- res$pc_loadings
  
  hotelling_table_df <- function(res,
                                 sample_col = "sample",
                                 class_col  = "class",
                                 digits_z   = 2L,
                                 digits_T2  = 2L,
                                 format     = FALSE) {
    ev <- res$extreme_values
    if (is.null(ev) || nrow(ev) == 0L) {
      return(ev[0, , drop = FALSE])
    }
    
    required_cols <- c(
      sample_col,
      class_col,
      "metabolite",
      "z_global",
      "abs_z_global",
      "z_class",
      "abs_z_class",
      "T2"
    )
    missing_cols <- setdiff(required_cols, names(ev))
    if (length(missing_cols) > 0L) {
      stop(
        "Missing columns in extreme_values: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    # Sort exactly like ui_outliers (but no head(top_n))
    ev_sorted <- ev[order(-ev$abs_z_global, -ev$abs_z_class, -ev$T2), , drop = FALSE]
    
    if (!format) {
      # Return numeric columns as-is (better for export / further analysis)
      out <- data.frame(
        Sample        = ev_sorted[[sample_col]],
        Class         = ev_sorted[[class_col]],
        Metabolite    = ev_sorted$metabolite,
        z_global      = ev_sorted$z_global,
        abs_z_global  = ev_sorted$abs_z_global,
        z_class       = ev_sorted$z_class,
        abs_z_class   = ev_sorted$abs_z_class,
        T2            = ev_sorted$T2,
        stringsAsFactors = FALSE
      )
    } else {
      # Return formatted character columns, matching the UI table styling
      out <- data.frame(
        Sample        = ev_sorted[[sample_col]],
        Class         = ev_sorted[[class_col]],
        Metabolite    = ev_sorted$metabolite,
        `Global z-score`   = formatC(ev_sorted$z_global,     format = "f", digits = digits_z),
        `|z| (global)`     = formatC(ev_sorted$abs_z_global, format = "f", digits = digits_z),
        `Class z-score`    = formatC(ev_sorted$z_class,      format = "f", digits = digits_z),
        `|z| (class)`      = formatC(ev_sorted$abs_z_class,  format = "f", digits = digits_z),
        `Mahalanobis^2`    = formatC(ev_sorted$T2,           format = "f", digits = digits_T2),
        stringsAsFactors = FALSE
      )
    }
    
    out
  }
  
  
  tab_numeric <- hotelling_table_df(res)
  
  shiny::withProgress(message = "Creating extreme_values_*today's_date*.xlsx...", value = 0, {
    s1 <- .add_sheet("Samples Outside Ellipse")
    txt1 <- paste(
      "Tab 1. This tab shows samples outside the Hotelling's T^2 95% ellipse.",
      "The ellipse is computed in the PC1-PC2 space using the non-QC samples in",
      "the signal drift corrected data."
    )
    openxlsx::writeData(wb,
                        s1,
                        x = txt1,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:20, rows = 1)
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
      x = outlier_samples,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 3, detail = "Saved: Samples Outside Ellipse")
    
    s2 <- .add_sheet("PC Loadings")
    txt2 <- paste(
      "Tab 2. This tab shows the loadings for PC1 and PC2 computed using the",
      "non-QC samples in the corrected data."
    )
    openxlsx::writeData(wb,
                        s2,
                        x = txt2,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s2, cols = 1:12, rows = 1)
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
      x = res$pc_loadings,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 3, detail = "Saved: PC Loadings")
    
    s3 <- .add_sheet("Potential Extreme Values")
    txt3 <- paste(
       "Tab 3. This tab shows samples outside the Hotelling's T^2 95% ellipse",
       "with AND have at least 1 potential extreme metabolite value meaning",
       "global AND class |z| is greater than 3 in the corrected data."
    )
    openxlsx::writeData(wb,
                        s3,
                        x = txt3,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s3, cols = 1:16, rows = 1)
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
      x = tab_numeric,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 3, detail = "Saved: Potential Extreme Values")
    
    if (!is.null(file)) {
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      return(normalizePath(file, winslash = "/"))
    }
  })
  return(wb)
}