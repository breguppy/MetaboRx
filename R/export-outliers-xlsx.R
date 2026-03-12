#' Export outliers excel file
#' @keywords internal
#' @noRd
export_outliers_xlsx <- function(p, d, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  
  wb <- openxlsx::createWorkbook()
  
  bold <- openxlsx::createStyle(textDecoration = "Bold")
  note <- openxlsx::createStyle(
    wrapText = TRUE,
    valign = "top",
    fgFill = "#f8cbad"
  )
  
  .add_sheet <- function(name) {
    nm <- gsub("[\\[\\]\\*\\?:/\\\\]", "_", name)
    nm <- substr(nm, 1L, 31L)
    openxlsx::addWorksheet(wb, nm)
    nm
  }
  
  .write_described_sheet <- function(wb,
                                     sheet_name,
                                     description,
                                     x,
                                     merge_cols = NULL,
                                     start_row_data = 3L,
                                     start_col_data = 1L,
                                     header_style = bold,
                                     note_style = note,
                                     note_row_height = 60) {
    s <- .add_sheet(sheet_name)
    
    openxlsx::writeData(
      wb,
      s,
      x = description,
      startCol = 1,
      startRow = 1
    )
    
    if (!is.null(merge_cols) && merge_cols >= 1L) {
      openxlsx::mergeCells(wb, s, cols = 1:merge_cols, rows = 1)
      openxlsx::addStyle(
        wb,
        s,
        style = note_style,
        rows = 1,
        cols = 1,
        gridExpand = TRUE
      )
      openxlsx::setRowHeights(wb, s, rows = 1, heights = note_row_height)
    }
    
    openxlsx::writeData(
      wb,
      s,
      x = x,
      startRow = start_row_data,
      startCol = start_col_data,
      headerStyle = header_style
    )
    
    invisible(s)
  }
  
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
    
    ev_sorted <- ev[order(-ev$abs_z_global, -ev$abs_z_class, -ev$T2), , drop = FALSE]
    
    if (!format) {
      out <- data.frame(
        Sample       = ev_sorted[[sample_col]],
        Class        = ev_sorted[[class_col]],
        Metabolite   = ev_sorted$metabolite,
        z_global     = ev_sorted$z_global,
        abs_z_global = ev_sorted$abs_z_global,
        z_class      = ev_sorted$z_class,
        abs_z_class  = ev_sorted$abs_z_class,
        T2           = ev_sorted$T2,
        stringsAsFactors = FALSE
      )
    } else {
      out <- data.frame(
        Sample             = ev_sorted[[sample_col]],
        Class              = ev_sorted[[class_col]],
        Metabolite         = ev_sorted$metabolite,
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
  
  df <- d$filtered_corrected$df_no_mv
  res <- detect_hotelling_nonqc_dual_z(df, p)
  
  outlier_samples <- unique(res$data$sample[res$data$is_outlier_sample])
  tab_numeric <- hotelling_table_df(res)
  
  shiny::withProgress(
    message = "Creating extreme_values_*today's_date*.xlsx...",
    value = 0,
    {
      .write_described_sheet(
        wb = wb,
        sheet_name = "Samples Outside Ellipse",
        description = paste(
          "Tab 1. This tab shows samples outside the Malahanois 95% ellipse.",
          "The ellipse is computed in the PC1-PC2 space using the non-QC samples in",
          "the signal drift corrected data."
        ),
        x = data.frame(Sample = outlier_samples, stringsAsFactors = FALSE),
        merge_cols = 12
      )
      shiny::incProgress(1 / 5, detail = "Saved: Samples Outside Ellipse")
      
      .write_described_sheet(
        wb = wb,
        sheet_name = "PC Loadings",
        description = paste(
          "Tab 2. This tab shows the loadings for PC1 and PC2 computed using the",
          "non-QC samples in the corrected data."
        ),
        x = res$pc_loadings,
        merge_cols = 12
      )
      shiny::incProgress(1 / 5, detail = "Saved: PC Loadings")
      
      .write_described_sheet(
        wb = wb,
        sheet_name = "Potential Extreme Values",
        description = paste(
          "Tab 3. This tab shows samples outside the Mahalanobis 95% ellipse",
          "that also have at least 1 potential extreme metabolite value, meaning",
          "global AND class |z| is greater than 3 in the corrected data."
        ),
        x = tab_numeric,
        merge_cols = 12
      )
      shiny::incProgress(1 / 5, detail = "Saved: Potential Extreme Values")
      
      global_export <- res$z_global
      global_export[setdiff(names(global_export), c("sample", "batch", "class", "order"))] <-
        lapply(
          global_export[setdiff(names(global_export), c("sample", "batch", "class", "order"))],
          function(x) round(x, 3)
        )
      
      .write_described_sheet(
        wb = wb,
        sheet_name = "Global Z-Scores",
        description = paste(
          "Tab 4. This tab shows pooled global z-scores for each retained metabolite",
          "in the corrected data. These z-scores are computed by centering and scaling",
          "all samples using the pooled non-QC samples only."
        ),
        x = global_export,
        merge_cols = min(12L, ncol(res$z_global))
      )
      shiny::incProgress(1 / 5, detail = "Saved: Global Z-Scores")
      
      class_export <- res$z_class
      class_export[setdiff(names(class_export), c("sample", "batch", "class", "order"))] <-
        lapply(
          class_export[setdiff(names(class_export), c("sample", "batch", "class", "order"))],
          function(x) round(x, 3)
        )
      
      .write_described_sheet(
        wb = wb,
        sheet_name = "Class Z-Scores",
        description = paste(
          "Tab 5. This tab shows class-based z-scores for each retained metabolite",
          "in the corrected data. For non-QC samples, z-scores are computed within",
          "their class using the class-specific mean and standard deviation.",
          "QC rows are expected to be NA on this sheet."
        ),
        x = class_export,
        merge_cols = min(12L, ncol(res$z_class))
      )
      shiny::incProgress(1 / 5, detail = "Saved: Class Z-Scores")
      
      if (!is.null(file)) {
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
        return(normalizePath(file, winslash = "/"))
      }
    }
  )
  
  wb
}