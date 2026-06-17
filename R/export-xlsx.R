#' exports corrected data as excel file
#'
#' @keywords internal
#' @noRd
export_xlsx <- function(p, d, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")

  wb <- openxlsx::createWorkbook()

  # make column names bold and descriptions with orange backgroud
  styles <- .xlsx_export_styles()
  bold <- styles$bold
  note <- styles$note

  .add_sheet <- function(name) {
    .xlsx_add_sheet(wb, name)
  }

  qc_rsd_flagged_metabolites <- function(df, rsd_cutoff) {
    if (!is.data.frame(df) ||
      length(rsd_cutoff) != 1L ||
      is.na(rsd_cutoff) ||
      !is.finite(rsd_cutoff)) {
      return(character(0))
    }

    rsd_df <- tryCatch(
      metabolite_rsd(df),
      error = function(e) NULL
    )

    if (is.null(rsd_df)) {
      return(character(0))
    }

    rsd_df$Metabolite[is.na(rsd_df$RSD_QC) | rsd_df$RSD_QC > rsd_cutoff]
  }

  add_eq_sheet <- identical(p$transform, "TRN")

  # base steps:
  # 1 Raw
  # 2 Settings
  # 3 Drift normalized
  # 4 Scaled/normalized
  # 5 Grouped organized
  # 6 Appendix1
  # +1 Fold change if control exists
  # +1 Appendix2 if extra metadata exists
  # +1 Equal weight sheet if TRN
  has_fold_change <- !isTRUE(p$no_control) && nzchar(p$control_class)

  meta_df <- d$cleaned$meta_df
  has_meta_appendix <- FALSE
  if (!is.null(meta_df)) {
    core_cols <- c("sample", "batch", "class", "order")
    extra_cols <- setdiff(names(meta_df), core_cols)
    has_meta_appendix <- length(extra_cols) > 0
  }

  N <- 6L +
    as.integer(has_fold_change) +
    as.integer(has_meta_appendix) +
    as.integer(add_eq_sheet)

  shiny::withProgress(
    message = "Creating corrected_data_*today's_date*.xlsx...",
    value = 0,
    {
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
      openxlsx::writeData(
        wb,
        s0,
        x = txt0,
        startCol = 1,
        startRow = 1
      )
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
          "Blank Threshold Multiplier",
          "Remove Blank-Threshold Failures?",
          "Sample/QC Average Difference Threshold",
          "Remove Sample/QC Average Failures?",
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
          sprintf(
            "%s%%",
            p$rsd_filter_threshold %||%
              p$rsd_cutoff %||%
              d$filtered_corrected$rsd_cutoff
          ),
          d$filtered$blank_threshold %||% NA,
          isTRUE(d$filtered$remove_blank_threshold_cols),
          sprintf("%s%%", d$filtered_corrected$percent_threshold %||% NA),
          isTRUE(p$remove_qc_average_pct_filter),
          p$transform,
          isTRUE(p$ex_ISTD),
          isTRUE(p$keep_corrected_qcs)
        )
      )
      txt1 <- paste(
        "Tab 1. This tab contains the correction settings along a list of metabolites",
        "that were eliminated throughout the process. The post-QC corrected data are",
        "shown on the next tab."
      )
      openxlsx::writeData(
        wb,
        s1,
        x = txt1,
        startCol = 1,
        startRow = 1
      )
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
      width_vec <- apply(
        settings,
        2,
        function(x) max(nchar(as.character(x)) + 2, na.rm = TRUE)
      )
      width_vec_header <- nchar(colnames(settings)) + 2
      max_width_vec <- pmax(width_vec, width_vec_header)
      openxlsx::setColWidths(wb, s1, cols = 1:2, widths = max_width_vec)

      if (isTRUE(p$remove_imputed)) {
        qc_rsd_removed_columns <- d$filtered_corrected$removed_metabolites_mv
        trn_withheld_columns <- d$transformed$withheld_cols_mv
      } else {
        qc_rsd_removed_columns <- d$filtered_corrected$removed_metabolites_no_mv
        trn_withheld_columns <- d$transformed$withheld_cols_no_mv
      }

      rsd_filter_disabled <- isTRUE(p$post_cor_filter) ||
        is.infinite(d$filtered_corrected$rsd_cutoff %||% Inf)
      qc_rsd_columns <- if (rsd_filter_disabled) {
        qc_rsd_flagged_metabolites(
          if (isTRUE(p$remove_imputed)) {
            d$filtered_corrected$df_mv
          } else {
            d$filtered_corrected$df_no_mv
          },
          p$rsd_filter_threshold %||%
            p$rsd_cutoff %||%
            d$filtered_corrected$rsd_cutoff
        )
      } else {
        qc_rsd_removed_columns
      }
      qc_rsd_header <- if (rsd_filter_disabled) {
        "QC-RSD Flagged Metabolites"
      } else {
        "QC-RSD Filtered Metabolites"
      }

      blank_threshold_columns <- if (isTRUE(d$filtered$remove_blank_threshold_cols)) {
        d$filtered$removed_blank_threshold_cols
      } else {
        d$filtered$blank_threshold_result$below_blank_threshold_ex_ISTD
      }
      blank_threshold_header <- if (isTRUE(d$filtered$remove_blank_threshold_cols)) {
        "Blank-Threshold Filtered Metabolites"
      } else {
        "Blank-Threshold Flagged Metabolites"
      }

      qc_average_columns <- if (isTRUE(p$remove_qc_average_pct_filter)) {
        d$filtered_corrected$removed_mets_pct_diff
      } else {
        d$filtered_corrected$flagged_mets
      }
      qc_average_header <- if (isTRUE(p$remove_qc_average_pct_filter)) {
        "Sample/QC Average Filtered Metabolites"
      } else {
        "Sample/QC Average Flagged Metabolites"
      }

      cur_col <- 4L
      add_settings_table <- function(cur_col, x, header = NULL) {
        if (is.null(x)) {
          return(cur_col)
        }

        if (is.data.frame(x)) {
          if (nrow(x) == 0L || ncol(x) == 0L) {
            return(cur_col)
          }
          df <- x
          if (!is.null(header)) {
            names(df) <- header
          }
        } else {
          vec <- as.character(stats::na.omit(unlist(x, use.names = FALSE)))
          if (length(vec) == 0L) {
            return(cur_col)
          }
          df <- stats::setNames(
            data.frame(vec, check.names = FALSE),
            header
          )
        }

        openxlsx::writeData(
          wb,
          s1,
          x = df,
          startRow = 3,
          startCol = cur_col,
          headerStyle = bold
        )
        width_vec <- vapply(
          seq_along(df),
          function(i) {
            max(nchar(c(names(df)[i], as.character(df[[i]]))) + 2, na.rm = TRUE)
          },
          numeric(1L)
        )
        openxlsx::setColWidths(
          wb,
          s1,
          cols = cur_col:(cur_col + ncol(df) - 1L),
          widths = width_vec
        )
        cur_col + ncol(df) + 1L
      }

      if (isTRUE(p$withhold_cols)) {
        cur_col <- add_settings_table(
          cur_col,
          d$cleaned$withheld_cols,
          "Columns Withheld From Correction"
        )
      }
      cur_col <- add_settings_table(
        cur_col,
        d$cleaned$all_missing_zero_qc_cols,
        "All Missing/Zero in QC Metabolites Removed"
      )
      cur_col <- add_settings_table(
        cur_col,
        d$cleaned$duplicate_col_names,
        "Duplicate Raw Column Names Repaired"
      )
      cur_col <- add_settings_table(
        cur_col,
        d$cleaned$duplicate_mets,
        c("Equal/Duplicate Metabolite Pairs Detected", "Paired Metabolite")
      )
      cur_col <- add_settings_table(
        cur_col,
        blank_threshold_columns,
        blank_threshold_header
      )
      cur_col <- add_settings_table(
        cur_col,
        d$filtered$mv_removed_cols,
        "Missing-Value Filtered Metabolites"
      )
      cur_col <- add_settings_table(
        cur_col,
        d$filtered$qc_missing_mets,
        "QC Missing After Missing-Value Filtering"
      )
      cur_col <- add_settings_table(
        cur_col,
        qc_average_columns,
        qc_average_header
      )
      cur_col <- add_settings_table(
        cur_col,
        qc_rsd_columns,
        qc_rsd_header
      )
      if (length(trn_withheld_columns) > 0) {
        cur_col <- add_settings_table(
          cur_col,
          trn_withheld_columns,
          "Excluded From Normalization"
        )
      }
      shiny::incProgress(1 / N, detail = "Saved: Correction Settings")

      # Sheet 2: Drift Normalized
      s2 <- .add_sheet("2. Drift Corrected")
      if (isTRUE(p$remove_imputed)) {
        df2 <- d$filtered_corrected$df_mv
      } else {
        df2 <- d$filtered_corrected$df_no_mv
      }
      if (!isTRUE(p$keep_corrected_qcs)) {
        df2 <- df2[df2$class != "QC", , drop = FALSE]
      }
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
        if (add_eq_sheet) {
          "to 1 for most metabolites. These data are equally weighted on the next tab before being normalized on the following tab."
        } else {
          "to 1 for most metabolites. These data are further normalized on the next tab."
        }
      )
      openxlsx::writeData(
        wb,
        s2,
        x = txt2,
        startCol = 1,
        startRow = 1
      )
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
      openxlsx::writeData(
        wb,
        s2,
        x = .round_metabolites_for_export(df2, .metabolite_cols_for_export(df2)),
        startRow = 3,
        headerStyle = bold
      )
      shiny::incProgress(1 / N, detail = "Saved: Drift Corrected")

      # Dynamic sheet numbering after Drift Normalized
      next_sheet_num <- 3L

      # Optional Sheet 3: Equal Weight (TRN only)
      if (add_eq_sheet) {
        s_eq <- .add_sheet(paste0(next_sheet_num, ". Equal Weight"))

        if (isTRUE(p$remove_imputed)) {
          df_eq <- d$transformed$equal_weight_df_mv
        } else {
          df_eq <- d$transformed$equal_weight_df_no_mv
        }

        if (!isTRUE(p$keep_corrected_qcs)) {
          df_eq <- df_eq[df_eq$class != "QC", , drop = FALSE]
        }

        txt_eq <- paste(
          "Tab", next_sheet_num, ". This tab shows metabolite level values after equal weighting.",
          "All metabolites included in TRN have been rescaled so that each metabolite has",
          "an average value of 1 across samples. This equally weights metabolites before",
          "total ratio normalization so that high-abundance metabolites do not dominate the",
          "normalization step."
        )

        openxlsx::writeData(
          wb,
          s_eq,
          x = txt_eq,
          startCol = 1,
          startRow = 1
        )
        openxlsx::mergeCells(wb, s_eq, cols = 1:16, rows = 1)
        openxlsx::addStyle(
          wb,
          s_eq,
          style = note,
          rows = 1,
          cols = 1,
          gridExpand = TRUE
        )
        openxlsx::setRowHeights(wb, s_eq, rows = 1, heights = 60)
        openxlsx::writeData(
          wb,
          s_eq,
          x = .round_metabolites_for_export(df_eq, .metabolite_cols_for_export(df_eq)),
          startRow = 3,
          headerStyle = bold
        )

        shiny::incProgress(1 / N, detail = "Saved: Equal Weight")
        next_sheet_num <- next_sheet_num + 1L
      }

      # Scaled or Normalized
      s_scaled <- .add_sheet(paste0(next_sheet_num, ". Samples Normalized"))
      df3 <- .samples_normalized_export_df(
        transformed = d$transformed,
        remove_imputed = p$remove_imputed,
        keep_corrected_qcs = p$keep_corrected_qcs
      )

      txt3 <- switch(p$transform,
        "TRN" = paste(
          "Tab", next_sheet_num, ". This tab shows metabolite level values",
          "ratiometrically normalized to total metabolite signal on a per",
          "sample basis. Before this normalization, metabolites were equally",
          "weighted so that each metabolite had an average value of 1000 across",
          "samples. This normalization is done by summing all individual",
          "post-QC corrected metabolite level values within a sample (total",
          "signal) and then dividing each individual metabolite level value",
          "within that sample by the total signal. These values are displayed",
          "on this tab after multiplying by the total number of metabolites",
          "present in the sample for easier visualization. Data remain in",
          "arbitrary units. Because arbitrary units for a given metabolite",
          "quantitatively scale across samples, levels of a given metabolite",
          "may be quantitatively compared across samples. Because unit scaling",
          "is different for each metabolite, different metabolites within a",
          "sample cannot be quantitatively compared. However, because",
          "differences in arbitrary unit scaling between samples cancel out by",
          "division, within-sample metabolite ratios can be quantitatively",
          "compared across samples."
        ),
        "ISTD_norm" = paste(
          "Tab", next_sheet_num, ". This tab shows the internal standard normalized",
          "metabolite level values. After correction, all metabolite levels are",
          "divided by the average of the internal standards within that sample."
        ),
        "none" = paste(
          "Tab", next_sheet_num, ". No scaling or normalization method has been applied to the data."
        )
      )

      openxlsx::writeData(
        wb,
        s_scaled,
        x = txt3,
        startCol = 1,
        startRow = 1
      )
      openxlsx::mergeCells(wb, s_scaled, cols = 1:22, rows = 1)
      openxlsx::addStyle(
        wb,
        s_scaled,
        style = note,
        rows = 1,
        cols = 1,
        gridExpand = TRUE
      )
      openxlsx::setRowHeights(wb, s_scaled, rows = 1, heights = 80)
      openxlsx::writeData(
        wb,
        s_scaled,
        x = .round_metabolites_for_export(df3, .metabolite_cols_for_export(df3)),
        startRow = 3,
        headerStyle = bold
      )
      shiny::incProgress(1 / N, detail = "Saved: Samples Normalized")

      next_sheet_num <- next_sheet_num + 1L

      # Grouped Data Organized
      s_grouped <- .add_sheet(paste0(next_sheet_num, ". Grouped Data Organized"))
      gdat <- group_stats(df3)
      openxlsx::writeData(
        wb,
        s_grouped,
        x = paste(
          "Tab", next_sheet_num, ". This tab shows post-scaled/normalized metabolite",
          "level values sorted by group, with group means, standard",
          "errors (SE), and coefficients of variation (CV) shown.",
          "Because the Metabolomics Core does not perform formal statistical",
          "analysis, these statistical analyses are shown for your",
          "convenience and quick appraisal of the data. For publication,",
          "data should be analyzed according to the standards of your",
          "field, including with the help of a statistician when needed."
        ),
        startCol = 1,
        startRow = 1
      )
      openxlsx::mergeCells(wb, s_grouped, cols = 1:12, rows = 1)
      openxlsx::addStyle(
        wb,
        s_grouped,
        style = note,
        rows = 1,
        cols = 1,
        gridExpand = TRUE
      )
      openxlsx::setRowHeights(wb, s_grouped, rows = 1, heights = 60)
      r <- 3
      for (nm in names(gdat$group_dfs)) {
        openxlsx::writeData(
          wb,
          s_grouped,
          x = .round_metabolites_for_export(
            gdat$group_dfs[[nm]],
            .metabolite_cols_for_export(gdat$group_dfs[[nm]])
          ),
          startRow = r,
          headerStyle = bold
        )
        r <- r + nrow(gdat$group_dfs[[nm]]) + 1
        openxlsx::writeData(
          wb,
          s_grouped,
          x = .round_metabolites_for_export(
            gdat$group_stats_dfs[[nm]],
            .metabolite_cols_for_export(gdat$group_stats_dfs[[nm]])
          ),
          startRow = r,
          startCol = 2,
          headerStyle = bold
        )
        r <- r + 6
      }
      shiny::incProgress(1 / N, detail = "Saved: Grouped Data Organized")

      next_sheet_num <- next_sheet_num + 1L

      # Grouped Data Fold Change
      if (has_fold_change) {
        s_fc <- .add_sheet(paste0(next_sheet_num, ". Grouped Data Fold Change"))
        ctrl_stats <- gdat$group_stats_dfs[[p$control_class]]
        fc_df <- fold_changes(df3, ctrl_stats[1, ])
        gfc <- group_stats(fc_df)
        r <- 1
        openxlsx::writeData(
          wb,
          s_fc,
          x = paste(
            "Tab", next_sheet_num, ". This tab shows post-scaled or normalized metabolite",
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
            "your field, including the help of a statistician when needed."
          ),
          startCol = 1,
          startRow = r
        )
        openxlsx::mergeCells(wb, s_fc, cols = 1:12, rows = r)
        openxlsx::addStyle(
          wb,
          s_fc,
          style = note,
          rows = 1,
          cols = 1,
          gridExpand = TRUE
        )
        openxlsx::setRowHeights(wb, s_fc, rows = r, heights = 60)
        r <- r + 2
        for (nm in names(gfc$group_dfs)) {
          openxlsx::writeData(
            wb,
            s_fc,
            x = .round_metabolites_for_export(
              gfc$group_dfs[[nm]],
              .metabolite_cols_for_export(gfc$group_dfs[[nm]])
            ),
            startRow = r,
            headerStyle = bold
          )
          r <- r + nrow(gfc$group_dfs[[nm]]) + 1
          openxlsx::writeData(
            wb,
            s_fc,
            x = .round_metabolites_for_export(
              gfc$group_stats_dfs[[nm]],
              .metabolite_cols_for_export(gfc$group_stats_dfs[[nm]])
            ),
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

      # Appendix1. MetaboAnalyst Ready
      s6 <- .add_sheet("Appendix1. MetaboAnalyst Ready")
      names(tf)[names(tf) == "sample"] <- "Sample Name"
      names(tf)[names(tf) == "class"] <- "Group"
      tf$batch <- NULL
      tf$order <- NULL
      openxlsx::writeData(
        wb,
        s6,
        x = .round_metabolites_for_export(tf, .metabolite_cols_for_export(tf))
      )
      shiny::incProgress(1 / N, detail = "Saved: MetaboAnalyst Ready")

      # Appendix2. MetaboAnalyst Meta tab (if extra metadata exists)
      if (has_meta_appendix) {
        core_cols <- c("sample", "batch", "class", "order")

        # Keep only rows present in Appendix1 (tf)
        meta_out <- meta_df[meta_df$sample %in% tf$`Sample Name`, , drop = FALSE]

        # enforce same ordering as tf
        meta_out <- meta_out[match(tf$`Sample Name`, meta_out$sample), , drop = FALSE]

        s7 <- .add_sheet("Appendix2. MetaboAnalyst Meta")

        # Rename to MetaboAnalyst convention
        names(meta_out)[names(meta_out) == "sample"] <- "Sample Name"
        names(meta_out)[names(meta_out) == "class"] <- "Group"

        # Remove unused columns
        meta_out$batch <- NULL
        meta_out$order <- NULL

        openxlsx::writeData(wb, s7, x = meta_out)
        shiny::incProgress(1 / N, detail = "Saved: MetaboAnalyst Meta")
      }

      if (!is.null(file)) {
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
        return(normalizePath(file, winslash = "/"))
      }
    }
  )

  return(wb)
}
