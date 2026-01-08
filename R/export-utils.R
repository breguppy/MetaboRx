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

#' Report text helpers for the HTML report
#'
#' @keywords internal
#' @noRd
report_text_withheld_columns <- function(p, d) {
  x <- htmltools::tagList(
    htmltools::tags$span(
      style = "font-weight:bold;",
      "The following columns are non-metabolite columns providing meta-information about the data:"
    ),
    htmltools::tags$ul(lapply(c(
      "sample = Identifies sample name",
      "batch = Identifies batch (large sample sets are separated into batches)",
      "class = Identifies sample type",
      "order = Identifies the order in which samples were injected into the instrument"
    ), htmltools::tags$li))
  )
  
  if (isTRUE(p$withhold_cols) && !is.null(p$n_withhold) && length(d$cleaned$withheld_cols) > 0) {
    x <- htmltools::tagList(
      x,
      htmltools::tags$br(),
      htmltools::tags$span(style = "font-weight:bold;", "The following columns were withheld from correction:"),
      htmltools::tags$ul(lapply(d$cleaned$withheld_cols, htmltools::tags$li))
    )
  }
  
  x
}

#' @keywords internal
#' @noRd
report_text_imputation <- function(p, d) {
  parts <- character(0)
  
  if (!identical(d$imputed$qc_str, "nothing to impute")) {
    parts <- c(parts, sprintf("Missing QC values are imputed with %s.", d$imputed$qc_str))
  } else {
    parts <- c(parts, "No missing QC values.")
  }
  
  if (!identical(d$imputed$sam_str, "nothing to impute")) {
    parts <- c(parts, sprintf("Missing sample values are imputed with %s.", d$imputed$sam_str))
  } else {
    parts <- c(parts, "No missing sample values.")
  }
  
  if ((!identical(d$imputed$qc_str, "nothing to impute") || !identical(d$imputed$sam_str, "nothing to impute")) &&
      isTRUE(p$remove_imputed)) {
    parts <- c(parts, "Imputed values are removed after correction.")
  }
  
  htmltools::tags$p(paste(parts, collapse = " "))
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
report_text_candidate_outliers <- function(p, d) {
  htmltools::tags$ul(lapply(c(
    "Possible extreme samples are detected by first grouping samples (QC vs non-QC or by class) and computing RSD.",
    "Metabolites with unstable QC RSD (greater than 30%) are not tested for extreme values.",
    "Robust z-scores are computed for each value within metabolite by median centering and scaling (MAD, IQR/1.349, SD, or 1) within each group.",
    "Candidate extreme values are non-QC sample-metabolite pairs with a z-score beyond the threshold of 4 for metabolites with stable QC RSD (<= 20%) or 5 for metabolites with borderline QC RSD (20% < QC RSD <= 30%).",
    "Each candidate is then confirmed with a test chosen by group size: Rosner/ESD for n > 25 (records a strength ratio), otherwise Dixon (if uniquely extreme and 3 <= n <= 30) or Grubbs (if extreme).",
    "Tied or ineligible cases can still be confirmed when the sample's squared Mahalanobis distance is flagged (md_only).",
    "Squared Mahalanobis distance is computed in the robust PC score space within each group.",
    "We retain PCs to reach at least 80% variance.",
    "Then a robust covariance is computed using MCD, OGK, shrinkage, or classical covariance depending on sample size and PCs retained.",
    "Confirmed candidates are possible extreme values; investigate before removing."
  ), htmltools::tags$li))
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
