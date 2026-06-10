#' Internal function for selecting metabolites for report scatter plots and
#' and list of metabolites that increase QC rsd after correction.
#'
#' @keywords internal
#' @noRd

.get_top_two <- function(p, d) {
  if (p$rsd_cal == "met") {
    top2 <- metabolite_rsd(d$filtered$df) |>
      dplyr::select(Metabolite, RSD_NonQC_before = RSD_NonQC) |>
      dplyr::inner_join(
        metabolite_rsd(d$filtered_corrected$df_no_mv) |>
          dplyr::select(Metabolite, RSD_NonQC_after = RSD_NonQC),
        by = "Metabolite"
      ) |>
      dplyr::mutate(decrease = RSD_NonQC_before - RSD_NonQC_after) |>
      dplyr::filter(is.finite(decrease)) |>
      dplyr::arrange(dplyr::desc(decrease)) |>
      dplyr::slice_head(n = 2) |>
      dplyr::pull(Metabolite)
  } else {
    top2 <- class_metabolite_rsd(d$filtered$df) |>
      dplyr::filter(class != "QC") |>
      dplyr::select(Metabolite, RSD_before = RSD) |>
      dplyr::inner_join(
        class_metabolite_rsd(d$filtered_corrected$df_no_mv) |>
          dplyr::filter(class != "QC") |>
          dplyr::select(Metabolite, RSD_after = RSD),
        by = "Metabolite"
      ) |>
      dplyr::mutate(decrease = RSD_before - RSD_after) |>
      dplyr::filter(is.finite(decrease)) |>
      dplyr::arrange(dplyr::desc(decrease)) |>
      dplyr::distinct(Metabolite, .keep_all = TRUE) |>
      dplyr::slice_head(n = 2) |>
      dplyr::pull(Metabolite)
  }
}

.increased_qc_rsd <- function(d) {
  increased_qc <- class_metabolite_rsd(d$filtered$df) |>
    dplyr::filter(class == "QC") |>
    dplyr::select(Metabolite, RSD_before = RSD) |>
    dplyr::inner_join(
      class_metabolite_rsd(d$filtered_corrected$df_no_mv) |>
        dplyr::filter(class == "QC") |>
        dplyr::select(Metabolite, RSD_after = RSD),
      by = "Metabolite"
    ) |>
    dplyr::filter(RSD_after > RSD_before) |>
    dplyr::arrange(dplyr::desc(RSD_after - RSD_before)) |>
    dplyr::pull(Metabolite)
}

#--------- Report text helper for popovers and HTML report

#' text for data structure and information requirements
#'
#' @keywords internal
#' @noRd
report_text_data_req <- function() {
  htmltools::tagList(
    htmltools::tags$p(
      htmltools::strong("1. Acceptable file formats: "),
      ".csv, .xls, or .xlsx"
    ),
    htmltools::tags$p(
      htmltools::tags$u("Note:"),
      " Raw data must be on the ",
      htmltools::tags$u("first sheet"),
      " of a .xls or .xlsx file."
    ),
    htmltools::tags$p(htmltools::strong("2. Required data formatting:")),
    htmltools::tags$ol(
      type = "a",
      htmltools::tags$li(htmltools::strong("Rows = samples"), " (can be in any order)"),
      htmltools::tags$li(
        htmltools::strong("Columns = non-metabolite columns and metabolites"),
        " (can be in any order)"
      ),
      htmltools::tags$li(htmltools::strong("Non-metabolite columns:")),
      htmltools::tags$ul(
        htmltools::tags$li(
          htmltools::tags$p(
            htmltools::strong("sample column (required): "),
            "Column that contains unique sample names."
          )
        ),
        htmltools::tags$li(
          htmltools::tags$p(
            htmltools::strong("batch column (optional): "),
            "Column that contains batch information if samples were run in batches."
          )
        ),
        htmltools::tags$li(
          htmltools::tags$p(
            htmltools::strong("class column (required): "),
            "Column that indicated the type of sample. Must contain QC samples labeled as 'NA', 'QC', 'Qc', or 'qc'. If data contains blank samples, label them as 'blank'."
          )
        ),
        htmltools::tags$li(
          htmltools::tags$p(
            htmltools::strong("injection order column (required): "),
            "Column that indicates injection order."
          )
        ),
        htmltools::tags$li(
          htmltools::tags$p(
            htmltools::strong("additional meta-information columns (optional): "),
            "Any remaining non-metabolite columns need to be specified."
          )
        )
      )
    ),
    htmltools::tags$p(
      htmltools::strong("3. Injection order must begin and end with QC samples: "),
      "Data (excluding blank samples) must begin and end with QC samples when sorted by injection order."
    ),
    htmltools::tags$p(
      htmltools::strong("4. Internal Standard Metabolites Must be Labeled: "),
      "internal standard metabolites must have a column name that begins with 'ISTD', or 'ITSD'."
    )
  )
}

#' Text for data cleaning and inspection step.
#' @keywords internal
#' @noRd
report_text_data_inspection <- function() {
  htmltools::tagList(
    htmltools::tags$p(
      "Before correction, the app performs a lightweight cleaning and runs a set of checks ",
      "confirming the dataset is formatted correctly. This step identifies common quality issues ",
      "by flagging invaild metabolite values, duplicate metabolites/columns, and possible contaminates."
    ),
    htmltools::tags$h4("Cleaning performed"),
    htmltools::tags$ul(
      htmltools::tags$li(
        htmltools::tags$strong("Standardizes metadata columns: "),
        "the selected columns are renamed to ",
        htmltools::tags$code("sample"), ", ",
        htmltools::tags$code("batch"), ", ",
        htmltools::tags$code("class"), ", and ",
        htmltools::tags$code("order"), "."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Removes non-metabolite columns withheld from correction: "),
        "if you chose to withhold additional columns."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Drops fully non-numeric metabolite columns: "),
        "columns that contain no numeric values."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Sorts rows by injection order: "),
        "all downstream checks reflect run order."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Standardizes QC labels: "),
        "QC samples are consistently recognized."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Detects and removes blank samples: "),
        "if present, blanks are excluded from correction."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Converts invalid metabolite entries to missing values (NA): "),
        htmltools::tags$ul(
          htmltools::tags$li("non-numeric entries (text, symbols, etc.) \u2192 NA"),
          htmltools::tags$li("exact zeros \u2192 NA")
        )
      )
    ),
    htmltools::tags$h4("Checks and summaries reported"),
    htmltools::tags$ul(
      htmltools::tags$li(
        htmltools::tags$strong("Dataset structure: "),
        "number of metabolite columns, number of samples and QC injections, number of batches, ",
        "and the list of non-QC classes."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Missingness: "),
        "total missing values (including values converted to NA during cleaning)."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("QC coverage by batch: "),
        "counts of QC injections in each batch."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Potential duplicate metabolites (informational): "),
        "flags metabolite column pairs that are equal or nearly equal across the same non-missing rows."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Blank-related flags (informational): "),
        "when blank samples are present, flags metabolites whose QC signal is low relative to blanks for review."
      )
    ),
    htmltools::tags$h4("How to use this section"),
    htmltools::tags$ul(
      htmltools::tags$li(
        "Use the structure counts to verify the correct columns were selected and parsed."
      ),
      htmltools::tags$li(
        "Use the warnings to spot issues that may affect correction (missingness, low QC coverage, or blank-related concerns)."
      ),
      htmltools::tags$li(
        "Flags in this section are primarily for transparency; features are not automatically removed unless explicitly stated."
      )
    )
  )
}

#' Withheld columns text  (if any columns are withheld)
#' @keywords internal
#' @noRd
report_text_withheld_columns <- function(p, d) {
  if (isTRUE(p$withhold_cols) && !is.null(p$n_withhold) && length(d$cleaned$withheld_cols) > 0) {
    htmltools::tagList(
      htmltools::tags$br(),
      htmltools::tags$span(
        style = "font-weight:bold;",
        "The following metabolite columns were withheld from correction:"
      ),
      htmltools::tags$ul(lapply(d$cleaned$withheld_cols, htmltools::tags$li))
    )
  } else {
    htmltools::tags$span(
      style = "font-weight:bold;",
      "All metabolite columns in the raw data were included in the correction."
    )
  }
}

#' Text describing missing value filtering.
#' @keywords internal
#' @noRd
report_text_mv_filter <- function() {
  htmltools::tagList(
    htmltools::tags$p(
      "Metabolites with missing value percentage above the ",
      "selected threshold for at least 1 sample class are removed from the ",
      " dataset. After filtering by missing value percentage, metabolites that ",
      "have at least 1 missing value for QC samples are displayed. Since missing",
      "values for QC samples is not common, further investigation is need to ",
      "determine if the value is truly not detected."
    ),
    htmltools::tags$p(
      "If a metabolite is missing for all samples in a single ",
      "class, a warning will appear stating the class and metabolite with all ",
      "missing values."
    ),
    htmltools::tags$h4("How to use this section"),
    htmltools::tags$ul(
      htmltools::tags$li(
        "Review the metabolite that are removed based on missing value percentage."
      ),
      htmltools::tags$li("Double check that QC missing values are truly missing in the raw data."),
      htmltools::tags$li(
        "Use the 'missing_value_counts.xlsx' to investigate patterns in missing values by viewing counts by sample, metabolite, batch, class, and class-metabolite."
      ),
    )
  )
}

#' Text explaining metabolite correlations and Pearson's r
#' @keywords internal
#' @noRd
report_text_correlations <- function() {
  htmltools::tagList(
    htmltools::tags$p(
      "To investigate linear relationships between metabolites, Pearson's r is computed for all pairs in the dataset. ",
      "A strong positive linear correlation (Pearson's r near 1) means that as one metabolite increases, the other metabolite consistently increases proportionally."
    ),
    htmltools::tags$p(
      "All pairwise correlations are computed, but only pairs with a strong positive linear correlations are displayed here. ",
      "To further investigate metabolite correlations view the Excel file '*_metabolite_correlations.xlsx'. ",
      "Take special note of any pair of metabolites what have strong correlation, but no biological explanation and investigate further if needed. "
    )
  )
}

#' Text describing imputation methods.
#' @keywords internal
#' @noRd
report_text_imputation <- function(p, d) {
  paragraphs <- character(0)

  if (!identical(d$imputed$qc_str, "nothing to impute")) {
    paragraphs <- c(
      paragraphs,
      sprintf("Missing QC values are imputed with %s.", d$imputed$qc_str)
    )
  } else {
    paragraphs <- c(
      paragraphs,
      "Since there are no missing QC values, no imputation is necessary."
    )
  }

  if (!identical(d$imputed$sam_str, "nothing to impute")) {
    paragraphs <- c(
      paragraphs,
      sprintf("Missing sample values are imputed with %s.", d$imputed$sam_str)
    )
  } else {
    paragraphs <- c(
      paragraphs,
      "Since there are no missing sample values, no imputation is necessary."
    )
  }

  if (
    (!identical(d$imputed$qc_str, "nothing to impute") ||
      !identical(d$imputed$sam_str, "nothing to impute")) &&
      isTRUE(p$remove_imputed)
  ) {
    paragraphs <- c(
      paragraphs,
      "Imputed values are removed after correction."
    )
  }

  htmltools::tagList(
    htmltools::tags$p(paste(paragraphs, collapse = " ")),
    htmltools::tags$strong("Imputation method descriptions:"),
    htmltools::tags$ul(lapply(
      c(
        "Metabolite median / mean: Across all samples.",
        "QC-metabolite median / mean: Across QC samples only.",
        "Class-metabolite median / mean: Across samples grouping by class.",
        "Minimum / half-minimum: Common for left-censored LC–MS data. Left-censored data occur when metabolite intensities fall below the instrument’s detection limit, so their exact values are unknown but known to be small.",
        "KNN: k-nearest neighbors imputation.",
        "Zero: Not recommended unless biologically justified."
      ),
      htmltools::tags$li
    ))
  )
}

#' Text describing correction methods
#' @keywords internal
#' @noRd
report_text_correction <- function(p, d) {
  htmltools::tags$p(sprintf(
    "Data was corrected using %s. For each metabolite, %s This model regresses peak areas in experimental samples, on an individual metabolite basis, against peak areas in pooled quality control samples.",
    d$corrected$str,
    d$corrected$parameters
  ))
}

#' @keywords internal
#' @noRd
report_text_correction_descriptions <- function() {
  htmltools::tagList(
    htmltools::tags$strong("QC-based signal drift correction methods:"),
    htmltools::tags$ul(
      htmltools::tags$li(
        htmltools::strong("Local constant regression (Nadaraya–Watson estimator): "),
        "Uses the weighted mean of nearby QC samples (vs injection order) to estimate drift and correct samples.",
        htmltools::tags$ul(
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Strengths: "),
              "Most stable and least complex method."
            )
          ),
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Weaknesses: "),
              "Can underfit real signal drift and leave residual trend."
            )
          )
        )
      ),
      htmltools::tags$li(
        htmltools::strong("Local linear regression: "),
        "Uses weighted local lines fit to QC samples (vs injection order) to estimate drift and correct samples.",
        htmltools::tags$ul(
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Strengths: "),
              "Stable and captures gradual increasing or decreasing drift trends."
            )
          ),
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Weaknesses: "),
              "Cannot capture strong curvature well and can chase QC noise when QCs are sparse."
            )
          )
        )
      ),
      htmltools::tags$li(
        htmltools::strong("Local polynomial regression (QC-RLSC / LOESS): "),
        "Uses weighted local polynomials (quadratic by default) fit to QC samples (vs injection order) to estimate drift and correct samples.",
        htmltools::tags$ul(
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Strengths: "),
              "Captures smooth nonlinear (curved) drift trends."
            )
          ),
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Weaknesses: "),
              "Can overfit with sparse QCs and typically performs poorly for abrupt, step-like drift."
            )
          )
        )
      ),
      htmltools::tags$li(
        htmltools::strong("Random forest (QC-RFSC): "),
        "Fits a random forest model on QC samples (QC intensity vs injection order) to estimate drift and correct samples.",
        htmltools::tags$ul(
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Strengths: "),
              "Flexible; can model irregular drift and abrupt changes."
            )
          ),
          htmltools::tags$li(
            htmltools::tags$p(
              htmltools::strong("Weaknesses: "),
              "Prefers many QCs (often ≥12–15). Highest overfitting risk because it can memorize QC noise. Less interpretable and does not enforce a smooth drift curve."
            )
          )
        )
      )
    ),
    htmltools::tags$p(
      htmltools::strong("General note: "),
      "All methods require QC samples that span the run (ideally at regular frequency). If QCs are sparse, clustered, or unstable, correction can be unreliable."
    )
  )
}

#' Text explaining how RSD is computed
#' @keywords internal
#' @noRd
report_text_rsd_cal <- function() {
  htmltools::tags$p(
    htmltools::strong("Relative Standard Deviation (RSD): "),
    "RSD = (SD / mean) × 100, ",
    "where SD is the standard deviation and mean is the average metabolite signal ",
    "across samples. Mean and standard deviation are computed for each metabolite ",
    "with missing values removed. RSD describes variability as a ",
    "percentage of the mean. Metabolites with QC RSD above the set threshold ",
    "are removed from the dataset and are listed here."
  )
}
#' Extreme value detection description.
#' @keywords internal
#' @noRd
report_text_ev_detection <- function() {
  htmltools::tagList(
    htmltools::tags$p(
      "This step flags potential extreme values using a 2D PCA / Mahalanobis ",
      "distance on non-QC samples. The squared Mahalanobis distance is computed ",
      "in the PC1–PC2 space using a PCA model fit on non-QC samples:"
    ),
    htmltools::tags$ol(
      htmltools::tags$li(
        htmltools::strong("Log and Scale Metabolites for non-QC samples: "),
        "Applies log2(x + 1), then standardizes using pooled non-QC samples only."
      ),
      htmltools::tags$li(
        htmltools::strong("PCA fit (non-QC only): "),
        "Fits PCA on pooled non-QC rows with complete metabolite data; uses PC1–PC2."
      ),
      htmltools::tags$li(
        htmltools::strong("Mahalanobis distance in PC space for all samples: "),
        "Projects all complete rows (QC + non-QC) into PC1–PC2 and computes a squared Mahalanobis distance."
      ),
      htmltools::tags$li(
        htmltools::strong("Ellipse cutoff: "),
        "Flags samples outside the (1 − α) ellipse using a χ² cutoff with df = 2 (default α = 0.05 → 95%)."
      ),
      htmltools::tags$li(
        htmltools::strong("Dual z-score rule (only for outlier samples): "),
        "Within samples outside the ellipse, flags metabolite values only when BOTH ",
        htmltools::strong("|global z| ≥ 3"),
        " (pooled non-QC scaling) AND ",
        htmltools::strong("|class z| ≥ 3"),
        " (within that non-QC class)."
      )
    ),
    htmltools::tags$hr(),
    htmltools::tags$p(
      htmltools::strong("Interpretation: "),
      "Red points are samples outside the ellipse. The table reports the specific metabolite values that also satisfy the dual z-score threshold."
    ),
    htmltools::tags$p(
      htmltools::tags$b("Caution: ", style = "color: red;"),
      "Candidate extreme values are displayed for the user's benefit. ",
      "Further investigation and justification is needed before categorizing an extreme value as an outlier and removing it."
    ),
    htmltools::tags$p(
      htmltools::strong("Note:"),
      "This is not an exhaustive outlier search. This extreme value detection sections aims to help identify unreasonable ",
      "metabolite values thats make a sample stand out from the others."
    )
  )
}

#' @keywords internal
#' @noRd
report_text_transformation <- function(p, d) {
  base <- htmltools::tags$p(d$transformed$str)

  if (length(d$transformed$withheld_cols) > 0) {
    return(htmltools::tagList(
      base,
      htmltools::tags$span(style = "font-weight:bold;", "The following columns are withheld from the transformation:"),
      htmltools::tags$ul(lapply(d$transformed$withheld_cols, htmltools::tags$li))
    ))
  }

  base
}

#' @keywords internal
#' @noRd
report_text_transform_methods <- function() {
  htmltools::tagList(
    htmltools::tags$p("Post-correction transformations/normalizations:"),
    htmltools::tags$ul(
      htmltools::tags$li(
        htmltools::strong("Internal Standard Normalization: "),
        "For each sample (row), compute the mean of internal standard columns (ISTD*/ITSD*), ",
        "then divide each non-internal standard metabolite by that mean."
      ),
      htmltools::tags$li(
        htmltools::strong("Total Ratio Normalization (TRN): "),
        "For each sample (row), compute the total signal as the sum across included metabolite columns, ",
        "then scale each included metabolite by its proportion of that total and multiply by the number ",
        "of non-missing metabolites in the sample (i.e., values become comparable across samples in arbitrary units)."
      ),
      htmltools::tags$li(
        htmltools::strong("None: "),
        "Leaves corrected metabolite values unchanged."
      ),
    ),
    htmltools::tags$hr(),
    htmltools::tags$p(
      htmltools::strong("Exclude internal standards checkbox: "),
      "When checked, ISTD/ITSD columns are excluded from the TRN total-signal calculation and are not transformed."
    ),
    htmltools::tags$p(
      htmltools::strong("Withhold from TRN: "),
      "Use this if a column should not contribute to the TRN total (e.g., TIC)."
    ),
    htmltools::tags$hr(),
    htmltools::tags$p(
      htmltools::tags$b("Caution: ", style = "color: red;"),
      "Internal Standard Normalization should not be used when only a single internal standard is measured.",
      "The single internal standard may not be representive of all metabolites measured in the samples.",
      "Total Ratio Normalization (TRN) relies on the assuption that total intensity should be the same across as samples",
      " and is sensitive to extreme values or dominate metabolites with high relative intensities."
    )
  )
}

#' How to interpret the metabolite scatter plots:
#' @keywords internal
#' @noRd
report_text_met_scatter <- function() {
  htmltools::tagList(
    htmltools::tags$p(htmltools::strong("How to read this plot")),
    htmltools::tags$p(
      "Dots are samples. In both panels, samples are organized in injection order on the x-axis with ",
      "QC samples colored blue and other samples colored yellow.",
      "Regardless of the correction method selected, the top panel of the figure shows ",
      "the selected metabolite's intensity values in the raw data. The bottom panel shows the metabolite's ",
      "scaled intensity values after correction. If mutiple batches are present in the dataset, the background ",
      "of the plot will alternate light gray and white to indicate different batches."
    ),
    htmltools::tags$p(
      htmltools::strong("Note: "),
      "the scale in the y-axes is different between the top and bottom plots."
    ),
    htmltools::tags$hr(),
    htmltools::tags$p(
      htmltools::strong("Correction method: Local polynomial fit (LOESS)", ),
      "For this correction method, the scatter plots will also show a smooth blue line. ",
      "The blue line is the local polynomial fit (LOESS) for QC samples and it summarizes QC signal drift over time. ",
      "The blue ribbon/shading around the line (if visible) shows uncertainty in that trend; narrow means ",
      "stable, wide means variable."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "after correction, the QC samples should have less variability and be centered around 1."
    ),
    htmltools::tags$hr(),
    htmltools::tags$p(
      htmltools::strong("Correction method: Random forest"),
      "For this correction method, the scatter plots will also show dashed and solid horizonal lines. ",
      "The tighter black dashed line is \u00B11 standard deviation (SD) around the QC mean. ",
      "The wider dark red solid line is \u00B12 SD around the QC mean. ",
      "These horizontal lines show how far QC values typically vary."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "After correction, QC points should be more stable (less drift over order) and more tightly ",
      "clustered within the SD bands compared with the raw panel."
    ),
  )
}

#' @keywords internal
#' @noRd
report_text_scatter_intro <- function(p, d) {
  parts <- c(
    "These plots show metabolites before and after signal drift correction before any transformation is applied.",
    "The two metabolites shown above have the largest decrease in sample variation.",
    "The change in variation was determined by calculating relative standard deviation (RSD) for each metabolite"
  )

  if (identical(p$rsd_cal, "class_met")) {
    parts <- c(parts, "grouping by sample class.")
  } else {
    parts <- c(parts, ".")
  }

  if (isTRUE(!p$post_cor_filter)) {
    parts <- c(parts, sprintf(
      "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above %s%%.",
      p$rsd_cutoff
    ))
  }

  htmltools::tags$p(paste(parts, collapse = " "))
}
#' RSD comparison plot explanation.
#' @keywords internal
#' @noRd
report_text_rsd_plots <- function() {
  htmltools::tagList(
    htmltools::tags$p(htmltools::strong("How to read this plot")),
    htmltools::tags$p(
      "This section compares relative standard deviation (RSD) in the corrected or transformed and corrected data ",
      "(depending on the setting selected under 'Compare raw data to') to the raw data. ",
      "RSD is computed by dividing the standard deviation of each metabolite by the mean of that metabolite and is expressed ",
      "as a percentage. RSD is computed for each metabolite for QC samples and non-QC samples separtely. RSD can also be computed for non-QC ",
      "samples grouping samples by class type (depending on the settings selected under 'Calculate RSD by')."
    ),
    htmltools::tags$hr(),
    htmltools::strong("Visualize changes in RSD by: Distrbution"),
    htmltools::tags$p(
      "The distributions of RSDs in non-QC samples is displayed in the left panel and the distribution of RSDs in ",
      "QC samples is displayed in the right panel. ",
      "The blue distribution is RSD in the raw data before any correction or transformations is applied. ",
      "The orange distrubution is RSD in the corrected or transformed and corrected data. "
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "after correction/transformation and correction, the orange distributions should be shifted to the left compared to the blue distributions. ",
      "The orange distribution for QC samples should be tall and skinny with the highest density near zero."
    ),
    htmltools::tags$hr(),
    htmltools::strong("Visualize changes in RSD by: Scatter Plot"),
    htmltools::tags$p(
      "In the scatter plot comparison the x-axis is RSD before correction/transformation and correction and the y-axis is RSD after. ",
      "RSDs for non-QC samples are displayed in the left panel and QC samples in the right panel. ",
      "Red dots indicate that RSD increased after correction/transformation and correction. ",
      "Gray dot indicate no change in RSD after correction/transformation and correction. ",
      "Green dots indicate a decrease in RSD after correction/transformation. ",
      "The percentages of increased, no change, and decreased RSDs are shown at the top of each panel."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "after correction/transofrmation and correction, the majority of RSDs should decrease for QC samples. ",
      "Non-QC sample RSDs may or may not decrease dramatically after correction/transformation and correction."
    ),
  )
}

#' @keywords internal
#' @noRd
report_text_rsd_intro <- function(p, d) {
  increased_qc <- .increased_qc_rsd(d)

  main <- sprintf(
    "In these plots, RSD is computed per %s after %s and compared to the raw data.",
    if (identical(p$rsd_cal, "class_met")) "sample class and metabolite" else "metabolite",
    if (identical(p$rsd_compare, "filtered_cor_data")) "correction" else "transformation and correction"
  )

  extra <- character(0)
  if (isTRUE(!p$post_cor_filter)) {
    extra <- c(extra, sprintf(
      "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above %s%%.",
      p$rsd_cutoff
    ))
  }

  x <- htmltools::tags$p(paste(c(main, extra), collapse = " "))

  if (length(increased_qc) > 0) {
    x <- htmltools::tagList(
      x,
      htmltools::tags$p(
        htmltools::HTML(sprintf(
          "The following metabolites increased QC RSD after correction:<br/>%s<br/>More investigation is needed to determine if these metabolites should be excluded.",
          paste(increased_qc, collapse = ", ")
        ))
      )
    )
  }

  x
}
#' Change in RSD table description.
#' @keywords internal
#' @noRd
report_text_rsd_table <- function() {
  htmltools::tagList(
    htmltools::strong("Performance Metric"),
    htmltools::tags$p("\u0394 RSD = RSD after correction \u2212 RSD before correction"),
    htmltools::tags$p(
      "The first table shows the median change in (\u0394) RSD for both QC samples and non-QC samples.",
      "\u0394 Metabolite RSD is computed for all non-QC samples and \u0394 Class-Metabolite RSD is computed by ",
      "grouping samples based on the 'class' column."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "After correction, RSD should decrease for both QC and non-QC samples. ",
      " In this situation, a more negative number is desirable for all \u0394 metrics."
    ),
    htmltools::tags$hr(),
    htmltools::strong("Post-correction Change"),
    htmltools::tags$p(
      "The second table shows the percentages of RSDs that increased or decreased after correction.",
      "Metabolite RSD is computed for all non-QC samples and Class-Metabolite RSD is computed by grouping samples ",
      "based on the 'class' column."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "After correction, RSD should decrease for both QC and non-QC samples. ",
      " Ideally, the percentage decreased should be much higher than the percentage increased."
    )
  )
}

#' Change in RSD table description.
#' @keywords internal
#' @noRd
report_text_rsd_tc_table <- function() {
  htmltools::tagList(
    htmltools::strong("Performance Metric"),
    htmltools::tags$p("\u0394 RSD = RSD after correction and transformation \u2212 RSD in raw data"),
    htmltools::tags$p(
      "The first table shows the median change in (\u0394) RSD for both QC samples and non-QC samples.",
      "\u0394 Metabolite RSD is computed for all non-QC samples and \u0394 Class-Metabolite RSD is computed by ",
      "grouping samples based on the 'class' column."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "After correction, RSD should decrease for both QC and non-QC samples. ",
      " In this situation, a more negative number is desirable for all \u0394 metrics."
    ),
    htmltools::tags$hr(),
    htmltools::strong("Post-transformation Change"),
    htmltools::tags$p(
      "The second table shows the percentages of RSDs that increased or decreased after correction and transformation.",
      "Metabolite RSD is computed for all non-QC samples and Class-Metabolite RSD is computed by grouping samples ",
      "based on the 'class' column."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "After correction and transformation, RSD should decrease for both QC and non-QC samples. ",
      " Ideally, the percentage decreased should be much higher than the percentage increased."
    )
  )
}

#' PCA plot description
#' @keywords internal
#' @noRd
report_text_pca_plots <- function() {
  htmltools::tagList(
    htmltools::tags$p(htmltools::strong("What is principal component analysis (PCA)?")),
    htmltools::tags$p(
      "PCA is a dimension reduction technique that projects the original data onto components that capture the maxium variance in the data. ",
      "Principal conponent 1 (PC1) represents the most variance in the data. After PC1, PC2 represents the most variance in the remaining ",
      "data."
    ),
    htmltools::tags$hr(),
    htmltools::strong("PCA score plots"),
    htmltools::tags$p(
      "The left panel is the 2D PC plot for the raw data and the right panel is the 2D PC plot for the corrected/transformed and corrected data. ",
      "The x-axis is PC1 and y-axis is PC2. The percentage in the parentheses on the axis labels is the variance explained for each conponent. ",
      "Dots in this figure represent samples."
    ),
    htmltools::tags$p(
      htmltools::strong("Goal: "),
      "after correction/transformation and correction, biological variation should dominate technical variation and signal drift should ",
      "not be visible in right panel. "
    ),
    htmltools::tags$ul(
      htmltools::tags$li(htmltools::strong("When coloring the plot by class: "), "QC samples should cluster together in the right panel."),
      htmltools::tags$li(htmltools::strong("When coloring the plot by batch or order: "), "there should be no distinct color patterns in the right panel if samples were run using a random injection ordering")
    ),
    htmltools::tags$hr(),
    htmltools::strong("PCA loading plots"),
    htmltools::tags$p(
      "The loading values show how much a metabolite contributes to that PC and the top 5 metabolites for each PC are shown below the PCA plot. ",
      "The magnitude of the loading corresponds to the metabolite's strength of correlation to that PC. ",
      "A metabolite with a large magnitude (close to 1 or -1) has a strong influence/contribution to that PC ",
      "and a metabolite with a small magnitude close to 0 has weak influence/contribution to that PC. ",
      "A positive loading (green) means that a high value in that metabolite corresponds to a high value in that PC. ",
      "A negative loading (red) means a high value in that metabolite corresponds to a low value in that PC."
    ),
  )
}
