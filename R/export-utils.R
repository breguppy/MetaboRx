#' Internal function for selecting metabolites for report scatter plots and
#' and list of metabolites that increase QC rsd after correction.
#'
#' @keywords internal
#' @noRd 

.get_top_two <- function(p, d) {
  if (p$rsd_cal == "met") {
    top2 <- metabolite_rsd(d$filtered$df)  %>%
      select(Metabolite, RSD_NonQC_before = RSD_NonQC) %>%
      inner_join(
        metabolite_rsd(d$filtered_corrected$df_no_mv) %>%
          select(Metabolite, RSD_NonQC_after = RSD_NonQC),
        by = "Metabolite"
      ) %>%
      mutate(decrease = RSD_NonQC_before - RSD_NonQC_after) %>%
      filter(is.finite(decrease)) %>%
      arrange(desc(decrease)) %>%
      slice_head(n = 2) %>%
      pull(Metabolite)
  } else {
    top2 <- class_metabolite_rsd(d$filtered$df) %>%
      filter(class != "QC") %>%
      select(Metabolite, RSD_before = RSD) %>%
      inner_join(
        class_metabolite_rsd(d$filtered_corrected$df_no_mv) %>%
          filter(class != "QC") %>%
          select(Metabolite, RSD_after = RSD),
        by = "Metabolite"
      ) %>%
      mutate(decrease = RSD_before - RSD_after) %>%
      filter(is.finite(decrease)) %>%
      arrange(desc(decrease)) %>%
      distinct(Metabolite, .keep_all = TRUE) %>%
      slice_head(n = 2) %>%
      pull(Metabolite)
  }
}

.increased_qc_rsd <- function(d) {
  increased_qc <- class_metabolite_rsd(d$filtered$df) %>%
    filter(class == "QC") %>%
    select(Metabolite, RSD_before = RSD) %>%
    inner_join(
      class_metabolite_rsd(d$filtered_corrected$df_no_mv) %>%
        filter(class == "QC") %>%
        select(Metabolite, RSD_after = RSD),
      by = "Metabolite"
    ) %>%
    filter(RSD_after > RSD_before) %>%
    arrange(desc(RSD_after - RSD_before)) %>%
    pull(Metabolite)
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

#' @keywords internal
#' @noRd
report_text_withheld_columns <- function(p, d) {
  if (isTRUE(p$withhold_cols) && !is.null(p$n_withhold) && length(d$cleaned$withheld_cols) > 0) {
    htmltools::tagList(
      htmltools::tags$br(),
      htmltools::tags$span(style = "font-weight:bold;", 
                           "The following metabolite columns were withheld from correction:"),
      htmltools::tags$ul(lapply(d$cleaned$withheld_cols, htmltools::tags$li))
    )
  } else {
    htmltools::tags$span(style = "font-weight:bold;", 
                         "All metabolite columns in the raw data were included in the correction.")
  }
}
#' @keywords internal
#' @noRd
report_test_mv_filter <- function() {
  htmltools::tagList(
    htmltools::tags$p("Metabolites with missing value percentage above the selected threshold are removed from the dataset.",
                "Use the 'missing_value_counts.xlsx' to investigate patterns in missing values by viewing sample, metabolite, batch, and class missing value counts."),
    htmltools::tags$p("Metabolites that remain in the dataset after filtering and have at least 1 missing value for QC samples are provided for diagnostic purposes.",
                "Since missing values for QC samples is not common, further investigation is need to determine if the value is truly not detected.")
  )
}

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

#' @keywords internal
#' @noRd
report_text_rsd_intro <- function(p, d) {
  increased_qc <- .increased_qc_rsd(d)
  
  main <- sprintf(
    "In these plots, green indicates RSD decreased after %s, red indicates RSD increased after %s, and gray indicates no change in RSD.",
    if (p$rsd_compare == "filtered_cor_data") "correction" else "correction and transformation",
    if (p$rsd_compare == "filtered_cor_data") "correction" else "correction and transformation"
  )
  
  extra <- character(0)
  if (identical(p$rsd_cal, "class_met")) {
    extra <- c(extra, "For these figures RSD is calculated for each metabolite grouping by sample class.")
  }
  
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

#' @keywords internal
#' @noRd
report_text_pca_intro <- function(p, d) {
  htmltools::tags$p(sprintf(
    "This PCA plot shows both the raw data and %s data colored by %s.",
    if (p$pca_compare == "filtered_cor_data") "corrected" else "corrected and transformed",
    p$color_col
  ))
}


#' @keywords internal
#' @noRd
report_text_hotelling_detection <- function(p, d) {
  if (identical(p$out_data, "filtered_cor_data")) {
    d_type <- "corrected data"
  } else {
    d_type <- "transformed and corrected data"
  }
  htmltools::tagList(
    
    htmltools::tags$p(
      paste("Candidate extreme samples are identified in", d_type),
      " using a two-dimensional principal component analysis (PCA) / Hotelling’s T² framework fit on non-QC samples. ",
      "Hotelling’s T² is computed as the squared Mahalanobis distance in PC1–PC2 space derived from a PCA model trained on non-QC samples only."
    ),
    
    htmltools::tags$ol(
      htmltools::tags$li(
        htmltools::tags$strong("Log and scale metabolites (non-QC only): "),
        "Metabolite intensities are transformed using log2(x + 1) and standardized using pooled non-QC samples."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("PCA fit (non-QC samples): "),
        "PCA is fit on non-QC samples with complete metabolite data, retaining PC1 and PC2."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Hotelling’s T² in PC space: "),
        "All complete samples (QC and non-QC) are projected into PC1–PC2 space, and a squared Mahalanobis distance (Hotelling’s T²) is computed."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Ellipse cutoff: "),
        "Samples falling outside the (1 − α) confidence ellipse are flagged using a χ² cutoff with 2 degrees of freedom (default α = 0.05, corresponding to 95%)."
      ),
      htmltools::tags$li(
        htmltools::tags$strong("Dual z-score rule for metabolite-level flags: "),
        "For samples outside the ellipse, individual metabolite values are flagged only if they satisfy both ",
        htmltools::tags$strong("|global z| ≥ 3"),
        " (scaled using pooled non-QC samples) and ",
        htmltools::tags$strong("|class z| ≥ 3"),
        " (scaled within the corresponding non-QC class)."
      )
    ),
    
    htmltools::tags$p(
      htmltools::tags$strong("Interpretation: "),
      "Samples outside the ellipse represent multivariate extremes. Reported candidate extreme values are metabolite measurements that additionally satisfy the dual z-score criterion."
    ),
    
    htmltools::tags$p(
      htmltools::tags$strong("Caution: "),
      "Candidate extreme values are provided for diagnostic purposes. ",
      "Additional biological, technical, or experimental context should be considered before classifying a value as an outlier and removing it."
    )
  )
}
