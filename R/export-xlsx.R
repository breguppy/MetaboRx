#' exports corrected data as excel file
#'
#' @keywords internal
#' @noRd
export_xlsx <- function(p, d, file = NULL) {
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
  
  # make progress bar
  N <- ifelse(isTRUE(p$no_control), 6, 5)
  shiny::withProgress(message = "Creating corrected_data_*today's_date*.xlsx...", value = 0, {
    # Sheet 0: Raw Data
    s0 <- .add_sheet("0. Raw Data")
    openxlsx::writeData(
      wb,
      s0,
      x = d$filtered$df,
      startRow = 3,
      headerStyle = bold
    )
    txt0 <- paste(
      "Tab 0. This tab contains raw peak areas for stable isotope-labeled",
      "(13C, deuterium) internal standards (pre-fix of 0-ISTD) and detected",
      "metabolites for both pooled quality control (QC) and experimental samples.",
      "Data within this tab are formatted for instrument drift correction using the",
      "multiply analyzed pooled QC sample. Batch, class, and order are identifiers",
      "for the correction methods, with batch and class set to 1 unless different",
      "batches or sample classes are being analyzed within the same run.",
      "The correction settings are shown on the next tab."
    )
    openxlsx::writeData(wb,
                        s0,
                        x = txt0,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s0, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s0,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s0, rows = 1, heights = 60)
    shiny::incProgress(1 / N, detail = "Saved: Raw Data")
    
    # Sheet 1: Correction settings
    s1 <- .add_sheet("1. Correction Settings")
    settings <- data.frame(
      Settings = c(
        "Sample Column Name",
        "Batch Column Name",
        "Class Column Name",
        "Injection Order Column Name",
        "Missing Value Threshold",
        "QC Imputation",
        "Sample Imputation",
        "Correction Method",
        "Remove Imputed After Correction?",
        "QC RSD% Threshold",
        "Scaling/Normalization",
        "Exclude ISTD in Scaling/Normalization",
        "Keep Corrected QCs"
      ),
      Values = c(
        p$sample_col,
        p$batch_col,
        p$class_col,
        p$order_col,
        sprintf("%s%%", d$filtered$mv_cutoff),
        d$imputed$qc_str,
        d$imputed$sam_str,
        d$corrected$str,
        isTRUE(p$remove_imputed),
        sprintf("%s%%", d$filtered_corrected$rsd_cutoff),
        p$transform,
        isTRUE(p$ex_ISTD),
        isTRUE(p$keep_corrected_qcs)
      ),
      stringsAsFactors = FALSE
    )
    txt1 <- paste(
      "Tab 1. This tab contains the correction settings along a list of metabolites",
      "that were eliminated throughout the process. The post-QC corrected data are",
      "shown on the next tab."
    )
    openxlsx::writeData(wb,
                        s1,
                        x = txt1,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:3, rows = 1)
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
      x = settings,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    width_vec <- apply(settings, 2, function(x)
      max(nchar(as.character(x)) + 2, na.rm = TRUE))
    width_vec_header <- nchar(colnames(settings)) + 2
    max_width_vec <- pmax(width_vec, width_vec_header)
    openxlsx::setColWidths(wb, s1, cols = 1:2, widths = max_width_vec)
    
    # Optional lists
    # columns withheld from correction: d$cleaned$withheld_cols
    # columns filtered by missing value: d$filtered$mv_removed_cols
    # columns filtered by QC RSD: d$filtered_corrected$removed_metabolites(_no)_mv
    # columns withheld from TRN
    if (isTRUE(p$remove_imputed)){
      qc_rsd_removed_columns <- d$filtered_corrected$removed_metabolites_mv
      trn_withheld_columns <- d$transformed$withheld_cols_mv
    } else {
      qc_rsd_removed_columns <- d$filtered_corrected$removed_metabolites_no_mv
      trn_withheld_columns <- d$transformed$withheld_cols_no_mv
    }
    cur_col <- 4
    add_list <- function(vec, header) {
      if (length(vec) == 0)
        return(invisible(NULL))
      df <- data.frame(setNames(list(vec), header))
      openxlsx::writeData(
        wb,
        s1,
        x = df,
        startRow = 3,
        startCol = cur_col,
        headerStyle = bold
      )
      width_vec <- apply(df, 2, function(x)
        max(nchar(as.character(x)) + 2, na.rm = TRUE))
      width_vec_header <- nchar(colnames(df)) + 2
      max_width_vec <- pmax(width_vec, width_vec_header)
      openxlsx::setColWidths(wb, s1, cols = cur_col, widths = max_width_vec)
      assign("cur_col", cur_col + 2, inherits = TRUE)
    }
    if (isTRUE(p$withhold_cols))
      add_list(d$cleaned$withheld_cols,
               "Columns Withheld From Correction")
    add_list(d$filtered$mv_removed_cols,
             "Missing-Value Filtered Metabolites")
    add_list(d$filtered$qc_missing_mets,
             "Metabolites with QC Missing Values")
    add_list(qc_rsd_removed_columns,
             "QC-RSD Filtered Metabolites")
    if (length(trn_withheld_columns) > 0)
      add_list(trn_withheld_columns,
               "Excluded From Normalization")
    shiny::incProgress(1 / N, detail = "Saved: Correction Settings")
    
    # Sheet 2: Drift Normalized
    s2 <- .add_sheet("2. Drift Normalized")
    if (isTRUE(p$remove_imputed)) {
      df2 <- d$filtered_corrected$df_mv
    } else {
      df2 <- d$filtered_corrected$df_no_mv
    }
    if (!isTRUE(p$keep_corrected_qcs))
      df2 <- df2[df2$class != "QC", , drop = FALSE]
    txt2 <- paste(
      "Tab 2. This tab shows instrument drift corrected values for metabolite levels",
      "in experimental samples. Data is corrected using",
      d$corrected$str,
      "For each metabolite, this method",
      d$corrected$parameters,
      "This model regresses",
      "peak areas in experimental samples, on an individual metabolite basis, against",
      "peak areas in pooled quality control samples. This corrects for normal instrument",
      "drift during the run. It produces relative metabolite level values in arbitrary",
      "units. For a given metabolite across the entire run, these values average close",
      "to 1 for most metabolites. These data are further normalized on the next tab."
    )
    openxlsx::writeData(wb,
                        s2,
                        x = txt2,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s2, cols = 1:16, rows = 1)
    openxlsx::addStyle(
      wb,
      s2,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s2, rows = 1, heights = 60)
    openxlsx::writeData(wb,
                        s2,
                        x = df2,
                        startRow = 3,
                        headerStyle = bold)
    shiny::incProgress(1 / N, detail = "Saved: Drift Normalized")
    
    # Sheet 3: Scaled or Normalized
    s3 <- .add_sheet("3. Scaled or Normalized")
    if (isTRUE(p$remove_imputed)) {
      df3 <- d$transformed$df_mv
    } else {
      df3 <- d$transformed$df_no_mv
    }
    
    keep <- setdiff(names(df3), trn_withheld_columns)
    df3 <- if (!isTRUE(p$keep_corrected_qcs))
      df3[df3$class != "QC", keep, drop = FALSE]
    else
      df3[, keep, drop = FALSE]
    txt3 <- switch(
      p$transform,
      "TRN"  = paste(
        "Tab 3. This tab shows metabolite level values",
        "ratiometrically normalized to total metabolite signal on a per",
        "sample basis. This normalization is done by summing all individual",
        "post-QC corrected metabolite level values within a sample (total",
        "signal) and then dividing each individual metabolite level value",
        "within that sample by the total signal. This normalization quantifies",
        "individual metabolite values across samples based on their proportion",
        "to total metabolite load, in arbitrary units, within each individual",
        "sample. These values are displayed on this tab after multiplying",
        "by the total number of metabolites present in the sample for easier",
        "visualization. Data remain in arbitrary units. Because arbitary",
        "units for a given metabolite quantitatively scale across samples,",
        "levels of a given metabolite may be quantitatively compared across",
        "samples. Because unit scaling is different for each metabolite,",
        "different metabolites within in a sample cannot be quantitatively",
        "compared. However, because differences in arbitrary unit scaling",
        "between samples cancel out by divsion, within-sample metabolite",
        "ratios can be quantitatively compared across samples."
      ),
      "ISTD_norm" = paste("Tab 3. This tab shows the internal standard normalized",
      "metabolite level values. After correction, all metabolite levels are",
      "divided by the average of the internal standards within that sample."
      ),
      "none" = "Tab 3. No scaling or normalization method has been applied to the data."
    )
    openxlsx::writeData(wb,
                        s3,
                        x = txt3,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s3, cols = 1:22, rows = 1)
    openxlsx::addStyle(
      wb,
      s3,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s3, rows = 1, heights = 80)
    openxlsx::writeData(wb,
                        s3,
                        x = df3,
                        startRow = 3,
                        headerStyle = bold)
    shiny::incProgress(1 / N, detail = "Saved: Scaled or Normalized")
    
    # Sheet 4: Grouped Data organized
    s4 <- .add_sheet("4. Grouped Data Organized")
    gdat <- group_stats(df3)
    openxlsx::writeData(
      wb,
      s4,
      x = paste(
        "Tab 4. This tab shows post-scaled/normalized metabolite",
        "level values sorted by group, with group means, standard",
        "erorrs (SE), and coefficients of variation (CV) shown.",
        "Because the Metabolomics Core does not perform formal statistical",
        "analysis, these statistical analyses are shown for your",
        "convenience and quick appraisal of the data. For publication,",
        "data should be analyzed according to the standards of your",
        "field, including with the help of a statistician when needed."
      ),
      startCol = 1,
      startRow = 1
    )
    openxlsx::mergeCells(wb, s4, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s4,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s4, rows = 1, heights = 60)
    r <- 3
    for (nm in names(gdat$group_dfs)) {
      openxlsx::writeData(
        wb,
        s4,
        x = gdat$group_dfs[[nm]],
        startRow = r,
        headerStyle = bold
      )
      r <- r + nrow(gdat$group_dfs[[nm]]) + 1
      openxlsx::writeData(
        wb,
        s4,
        x = gdat$group_stats_dfs[[nm]],
        startRow = r,
        startCol = 2,
        headerStyle = bold
      )
      r <- r + 6
    }
    shiny::incProgress(1 / N, detail = "Saved: Grouped Data Organized")
    
    # Sheet 5: Grouped Data Fold Change
    if (!isTRUE(p$no_control) && nzchar(p$control_class)) {
      s5 <- .add_sheet("5. Grouped Data Fold Change")
      ctrl_stats <- gdat$group_stats_dfs[[p$control_class]]
      fc_df <- fold_changes(df3, ctrl_stats[1, ])
      gfc   <- group_stats(fc_df)
      r <- 1
      openxlsx::writeData(
        wb,
        s5,
        x = paste(
          "Tab 5. This tab shows post-scaled or normalized metabolite",
          "level values expressed in terms of fold change relative",
          "to the",
          p$control_class,
          "group mean. These values have",
          "been sorted by group, with group means, standard errors",
          "(SE), and coefficients of variation (CV) shown. Because",
          "the Metabolomics Core does not perform formal statistical",
          "analysis, these statistical analyses are shown for your",
          "convenience and quick appraisal of the data. For publication,",
          "data should be analyzed according to the standards of",
          "your field, including with the help of a statistician when needed."
        ),
        startCol = 1,
        startRow = r
      )
      openxlsx::mergeCells(wb, s5, cols = 1:12, rows = r)
      openxlsx::addStyle(
        wb,
        s5,
        style = note,
        rows = 1,
        cols = 1,
        gridExpand = TRUE
      )
      openxlsx::setRowHeights(wb, s5, rows = r, heights = 60)
      r <- r + 2
      for (nm in names(gfc$group_dfs)) {
        openxlsx::writeData(
          wb,
          s5,
          x = gfc$group_dfs[[nm]],
          startRow = r,
          headerStyle = bold
        )
        r <- r + nrow(gfc$group_dfs[[nm]]) + 1
        openxlsx::writeData(
          wb,
          s5,
          x = gfc$group_stats_dfs[[nm]],
          startRow = r,
          startCol = 2,
          headerStyle = bold
        )
        r <- r + 6
      }
      tf <- fc_df
      shiny::incProgress(1 / N, detail = "Saved: Grouped Data Fold Change")
    } else {
      tf <- df3
    }

    # Appendix1. Metaboanalyst Ready
    s6 <- .add_sheet("Appendix1. MetaboAnalyst Ready")
    names(tf)[names(tf) == "sample"] <- "Sample Name"
    names(tf)[names(tf) == "class"]  <- "Group"
    tf$batch <- NULL
    tf$order <- NULL
    openxlsx::writeData(wb, s6, x = tf)
    
    shiny::incProgress(1 / N, detail = "Saved: Metaboanalyst Ready")
    
    # Appendix2. MetaboAnalyst Meta tab (if extra metadata exists)
    meta_df <- d$cleaned$meta_df
    
    if (!is.null(meta_df)) {
      core_cols <- c("sample", "batch", "class", "order")
      extra_cols <- setdiff(names(meta_df), core_cols)
      
      if (length(extra_cols) > 0) {
        # Keep only rows present in Appendix1 (tf)
        meta_out <- meta_df[meta_df$sample %in% tf$`Sample Name`, , drop = FALSE]
        
        # Optional: enforce same ordering as tf
        meta_out <- meta_out[match(tf$`Sample Name`, meta_out$sample), , drop = FALSE]
        
        s7 <- .add_sheet("Appendix2. MetaboAnalyst Meta")
        
        # Rename to MetaboAnalyst convention
        names(meta_out)[names(meta_out) == "sample"] <- "Sample Name"
        names(meta_out)[names(meta_out) == "class"]  <- "Group"
        
        # Remove unused columns
        meta_out$batch <- NULL
        meta_out$order <- NULL
        
        openxlsx::writeData(wb, s7, x = meta_out)
        
        shiny::incProgress(1 / N, detail = "Saved: MetaboAnalyst Meta")
      }
    }
    
    if (!is.null(file)) {
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      return(normalizePath(file, winslash = "/"))
    }
  })
  return(wb)
}
