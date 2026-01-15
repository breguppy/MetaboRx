#' exports missing value report as excel file
#'
#' @keywords internal
#' @noRd
export_mv_xlsx <- function(p, d, file = NULL) {
  cleaned_df <- d$cleaned$df
  has_all_missing <- !is.null(d$filtered$class_metab_all_missing) &&
    is.data.frame(d$filtered$class_metab_all_missing) &&
    nrow(d$filtered$class_metab_all_missing) > 0L
  
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  
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
  
  shiny::withProgress(message = "Creating missing value summary...", value = 0, {
    meta_cols <- c("sample", "batch", "class", "order")
    metab_cols <- setdiff(names(cleaned_df), meta_cols)
    n_metabs <- length(metab_cols)
    
    # count missing values by metabolite (num missing in samples and missing in group and as percentages of data)
    count_missing_by_metabolite <- function(df, metab_cols) {
      n_samples <- nrow(df)
      
      tibble::tibble(
        metabolite    = metab_cols,
        missing_count = vapply(df[metab_cols], function(col) sum(is.na(col)), integer(1)),
        missing_pct   = (missing_count / n_samples) * 100
      ) |>
        dplyr::arrange(dplyr::desc(missing_pct))
    }
    
    sample_df <- cleaned_df[cleaned_df$class != "QC", ]
    qc_df <- cleaned_df[cleaned_df$class == "QC", ]
    
    sample_met_mv <- count_missing_by_metabolite(sample_df, metab_cols)
    qc_met_mv <- count_missing_by_metabolite(qc_df, metab_cols)
    
    sample_renamed <- sample_met_mv |>
      dplyr::rename(
        sample_missing_count = missing_count,
        sample_missing_pct   = missing_pct
      )
    
    qc_renamed <- qc_met_mv |>
      dplyr::rename(
        qc_missing_count = missing_count,
        qc_missing_pct   = missing_pct
      )
    
    metab_mv <- sample_renamed |>
      dplyr::inner_join(qc_renamed, by = "metabolite") |>
      dplyr::filter(!(sample_missing_pct == 0 & qc_missing_pct == 0)) |>
      dplyr::arrange(dplyr::desc(sample_missing_pct))
    
    # save them to workbook:
    # sheet 1 Metabolite
    s1 <- .add_sheet("Metabolite")
    txt1 <- paste(
      "Tab Metabolite. Missing value counts (missing_count) and percentages (missing_pct) per metabolite are listed here for samples and QC samples.",
      "If a metabolite is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb,
                        s1,
                        x = txt1,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:6, rows = 1)
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
      x = metab_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 4, detail = "Saved: missing values by metabolite")
  
    # count missing values by sample (num missing and percentage)
    sample_mv <- cleaned_df |>
      dplyr::mutate(
        missing_count = dplyr::across(dplyr::all_of(metab_cols), is.na) |>
          rowSums(),
        missing_pct   = (missing_count / n_metabs) * 100
      ) |>
      dplyr::select(sample, missing_count, missing_pct) |>
      dplyr::filter(missing_pct > 0) |>
      dplyr::arrange(dplyr::desc(missing_pct))
    
    # sheet 2 sample
    s2 <- .add_sheet("Sample")
    txt2 <- paste(
      "Tab Sample. Missing value counts (missing_count) and percentages (missing_pct) per sample are listed here.",
      "If a sample is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb,
                        s2,
                        x = txt2,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s2, cols = 1:6, rows = 1)
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
      x = sample_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 4, detail = "Saved: missing values by sample")
  
    # count missing values by class (num missing and as percentage)
    class_mv <- cleaned_df |>
      dplyr::group_by(class) |>
      dplyr::summarise(
        n_samples      = dplyr::n(),
        missing_count  = dplyr::across(
          dplyr::all_of(metab_cols),
          ~ sum(is.na(.x)),
          .unpack = FALSE
        ) |>
          unlist() |>
          sum(),
        total_values   = n_samples * n_metabs,
        missing_pct    = (missing_count / total_values) * 100,
        .groups = "drop"
      ) |>
      dplyr::select(class, missing_count, missing_pct) |>
      dplyr::filter(missing_pct > 0) |>
      dplyr::arrange(dplyr::desc(missing_pct))
    
    # sheet 3 class
    s3 <- .add_sheet("Class")
    txt3 <- paste(
      "Tab Class. Missing value counts (missing_count) and percentages (missing_pct) per sample class are listed here.",
      "If a sample class is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb,
                        s3,
                        x = txt3,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s3, cols = 1:6, rows = 1)
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
      x = class_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 4, detail = "Saved: missing values by class")
  
    # count missing by batch (num missing and as a percentage)
    batch_mv <- cleaned_df |>
      dplyr::group_by(batch) |>
      dplyr::summarise(
        n_samples      = dplyr::n(),
        missing_count  = dplyr::across(
          dplyr::all_of(metab_cols),
          ~ sum(is.na(.x)),
          .unpack = FALSE
        ) |>
          unlist() |>
          sum(),
        total_values   = n_samples * n_metabs,
        missing_pct    = (missing_count / total_values) * 100,
        .groups = "drop"
      ) |>
      dplyr::select(batch, missing_count, missing_pct) |>
      dplyr::filter(missing_pct > 0) |>
      dplyr::arrange(dplyr::desc(missing_pct))
    
    # sheet 4 batch
    s4 <- .add_sheet("Batch")
    txt4 <- paste(
      "Tab Batch. Missing value counts (missing_count) and percentages (missing_pct) per batch are listed here.",
      "If a batch is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb,
                        s4,
                        x = txt4,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s4, cols = 1:6, rows = 1)
    openxlsx::addStyle(
      wb,
      s4,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s4, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s4,
      x = batch_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 4, detail = "Saved: missing values by batch")
    
    if (has_all_missing) {
      s5 <- .add_sheet("All Missing for Class")
      txt5 <- paste(
        "Tab All Missing for Class. The following class-metabolite pairs have",
        "all values missing."
      )
      openxlsx::writeData(wb,
                          s5,
                          x = txt5,
                          startCol = 1,
                          startRow = 1)
      openxlsx::mergeCells(wb, s5, cols = 1:6, rows = 1)
      openxlsx::addStyle(
        wb,
        s5,
        style = note,
        rows = 1,
        cols = 1, 
        gridExpand = TRUE
      )
      openxlsx::writeData(
        wb,
        s5,
        x = d$filtered$class_metab_all_missing[, c("class", "metabolite"), drop = FALSE],
        startRow = 3,
        startCol = 1,
        headerStyle = bold
      )
    }
    
  })
  return(wb)
}