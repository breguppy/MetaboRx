#' exports missing value report as excel file
#'
#' @keywords internal
#' @noRd
export_mv_xlsx <- function(p, d, file = NULL) {
  cleaned_df <- d$cleaned$df
  
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  
  # make column names bold and descriptions with orange background
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
  
  shiny::withProgress(message = "Creating missing value summary...", value = 0, {
    meta_cols <- c("sample", "batch", "class", "order")
    metab_cols <- setdiff(names(cleaned_df), meta_cols)
    n_metabs <- length(metab_cols)
    
    # count missing values by metabolite (samples + QC)
    count_missing_by_metabolite <- function(df, metab_cols) {
      n_samples <- nrow(df)
      
      tibble::tibble(
        metabolite = metab_cols,
        missing_count = vapply(
          df[metab_cols],
          function(col) sum(is.na(col) | col <= 0),
          integer(1)
        ),
        missing_pct = if (n_samples == 0L) NA_real_ else (missing_count / n_samples) * 100
      ) |>
        dplyr::arrange(dplyr::desc(missing_pct))
    }
    
    sample_df <- cleaned_df[cleaned_df$class != "QC", , drop = FALSE]
    qc_df <- cleaned_df[cleaned_df$class == "QC", , drop = FALSE]
    
    sample_met_mv <- count_missing_by_metabolite(sample_df, metab_cols)
    qc_met_mv <- count_missing_by_metabolite(qc_df, metab_cols)
    
    sample_renamed <- sample_met_mv |>
      dplyr::rename(
        sample_missing_count = missing_count,
        sample_missing_pct = missing_pct
      )
    
    qc_renamed <- qc_met_mv |>
      dplyr::rename(
        qc_missing_count = missing_count,
        qc_missing_pct = missing_pct
      )
    
    metab_mv <- sample_renamed |>
      dplyr::inner_join(qc_renamed, by = "metabolite") |>
      dplyr::filter(!(sample_missing_pct == 0 & qc_missing_pct == 0)) |>
      dplyr::arrange(dplyr::desc(sample_missing_pct))
    
    # sheet 1 Metabolite
    s1 <- .add_sheet("Metabolite")
    txt1 <- paste(
      "Tab Metabolite. Missing value counts (missing_count) and percentages (missing_pct)",
      "per metabolite are listed here for samples and QC samples. Missing values are",
      "defined as NA or <= 0. If a metabolite is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb, s1, x = txt1, startCol = 1, startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:6, rows = 1)
    openxlsx::addStyle(wb, s1, style = note, rows = 1, cols = 1, gridExpand = TRUE)
    openxlsx::setRowHeights(wb, s1, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s1,
      x = metab_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by metabolite")
    
    # count missing values by sample
    sample_mv <- cleaned_df |>
      dplyr::mutate(
        missing_count = rowSums(dplyr::across(dplyr::all_of(metab_cols), ~ is.na(.x) | .x <= 0)),
        missing_pct = if (n_metabs == 0L) NA_real_ else (missing_count / n_metabs) * 100
      ) |>
      dplyr::select(sample, missing_count, missing_pct) |>
      dplyr::filter(missing_pct > 0) |>
      dplyr::arrange(dplyr::desc(missing_pct))
    
    # sheet 2 Sample
    s2 <- .add_sheet("Sample")
    txt2 <- paste(
      "Tab Sample. Missing value counts (missing_count) and percentages (missing_pct)",
      "per sample are listed here. Missing values are defined as NA or <= 0.",
      "If a sample is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb, s2, x = txt2, startCol = 1, startRow = 1)
    openxlsx::mergeCells(wb, s2, cols = 1:6, rows = 1)
    openxlsx::addStyle(wb, s2, style = note, rows = 1, cols = 1, gridExpand = TRUE)
    openxlsx::setRowHeights(wb, s2, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s2,
      x = sample_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by sample")
    
    # count missing values by class
    class_mv <- cleaned_df |>
      dplyr::group_by(class) |>
      dplyr::summarise(
        n_samples = dplyr::n(),
        missing_count = sum(unlist(dplyr::across(dplyr::all_of(metab_cols), ~ sum(is.na(.x) | .x <= 0)))),
        total_values = n_samples * n_metabs,
        missing_pct = ifelse(total_values == 0, NA_real_, (missing_count / total_values) * 100),
        .groups = "drop"
      ) |>
      dplyr::select(class, missing_count, missing_pct) |>
      dplyr::filter(missing_pct > 0) |>
      dplyr::arrange(dplyr::desc(missing_pct))
    
    # sheet 3 Class
    s3 <- .add_sheet("Class")
    txt3 <- paste(
      "Tab Class. Missing value counts (missing_count) and percentages (missing_pct)",
      "per sample class are listed here. Missing values are defined as NA or <= 0.",
      "If a sample class is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb, s3, x = txt3, startCol = 1, startRow = 1)
    openxlsx::mergeCells(wb, s3, cols = 1:6, rows = 1)
    openxlsx::addStyle(wb, s3, style = note, rows = 1, cols = 1, gridExpand = TRUE)
    openxlsx::setRowHeights(wb, s3, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s3,
      x = class_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by class")
    
    # count missing by batch
    batch_mv <- cleaned_df |>
      dplyr::group_by(batch) |>
      dplyr::summarise(
        n_samples = dplyr::n(),
        missing_count = sum(unlist(dplyr::across(dplyr::all_of(metab_cols), ~ sum(is.na(.x) | .x <= 0)))),
        total_values = n_samples * n_metabs,
        missing_pct = ifelse(total_values == 0, NA_real_, (missing_count / total_values) * 100),
        .groups = "drop"
      ) |>
      dplyr::select(batch, missing_count, missing_pct) |>
      dplyr::filter(missing_pct > 0) |>
      dplyr::arrange(dplyr::desc(missing_pct))
    
    # sheet 4 Batch
    s4 <- .add_sheet("Batch")
    txt4 <- paste(
      "Tab Batch. Missing value counts (missing_count) and percentages (missing_pct)",
      "per batch are listed here. Missing values are defined as NA or <= 0.",
      "If a batch is not listed here, it did not have any missing values."
    )
    openxlsx::writeData(wb, s4, x = txt4, startCol = 1, startRow = 1)
    openxlsx::mergeCells(wb, s4, cols = 1:6, rows = 1)
    openxlsx::addStyle(wb, s4, style = note, rows = 1, cols = 1, gridExpand = TRUE)
    openxlsx::setRowHeights(wb, s4, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s4,
      x = batch_mv,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by batch")
    
    # class-metabolite summary
    class_metab_missing_summary <- if (length(metab_cols) == 0L) {
      data.frame(
        class = character(0),
        metabolite = character(0),
        n_rows_in_class = integer(0),
        missing_count = integer(0),
        missing_pct = numeric(0),
        all_missing = logical(0)
      )
    } else {
      do.call(
        rbind,
        lapply(sort(unique(cleaned_df$class)), function(cl) {
          idx <- which(cleaned_df$class == cl)
          n_in_class <- length(idx)
          
          do.call(
            rbind,
            lapply(metab_cols, function(met) {
              x <- cleaned_df[idx, met]
              miss <- is.na(x) | x <= 0
              
              data.frame(
                class = cl,
                metabolite = met,
                n_rows_in_class = n_in_class,
                missing_count = sum(miss),
                missing_pct = if (n_in_class == 0L) NA_real_ else (sum(miss) / n_in_class) * 100,
                all_missing = if (n_in_class == 0L) TRUE else all(miss),
                row.names = NULL
              )
            })
          )
        })
      )
    }
    
    if (nrow(class_metab_missing_summary) > 0L) {
      s5 <- .add_sheet("Class-Met Missing")
      txt5 <- paste(
        "Tab Class-Met Missing. Missing value counts and percentages are listed",
        "for each class-metabolite pair. Missing values are defined as NA or <= 0.",
        "The column all_missing indicates whether all values for that metabolite",
        "within the class are missing."
      )
      openxlsx::writeData(wb, s5, x = txt5, startCol = 1, startRow = 1)
      openxlsx::mergeCells(wb, s5, cols = 1:6, rows = 1)
      openxlsx::addStyle(wb, s5, style = note, rows = 1, cols = 1, gridExpand = TRUE)
      openxlsx::setRowHeights(wb, s5, rows = 1, heights = 75)
      
      class_metab_missing_export <- class_metab_missing_summary |>
        dplyr::arrange(dplyr::desc(all_missing), dplyr::desc(missing_pct), class, metabolite)
      
      openxlsx::writeData(
        wb,
        s5,
        x = class_metab_missing_export,
        startRow = 3,
        startCol = 1,
        headerStyle = bold
      )
    }
    
    shiny::incProgress(1 / 5, detail = "Saved: missing values by class and metabolite")
  })
  
  return(wb)
}