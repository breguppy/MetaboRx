#' Writes HTML of correction report then writes PDF of correction report
#'
#' @keywords internal
#' @noRd
render_report <- function(p,
                          d,
                          out_dir,
                          template = system.file("app", "report_templates", "report.Rmd", package = "QCcorrection")) {
  .require_pkg("rmarkdown", "render reports")
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)
  env <- new.env(parent = baseenv())
  shiny::withProgress(message = "Creating quality_report.html...", value = 0, {
    # build plots
    met_candidates <- .get_top_two(p, d)
    increased_qc <- .increased_qc_rsd(d)
    met1_plot <- make_met_scatter(d, p, met_candidates[1])
    met2_plot <- make_met_scatter(d, p, met_candidates[2])
    rsd_plot  <- make_rsd_plot(p, d)
    pca_plot  <- make_pca_plot(p, d)
    pca_loading_plot <- make_pca_loading_plot(p, d)
    shiny::incProgress(1 / 3, detail = "Saved: plots for report")
    
    # plain strings only; no shiny::tagList here
    descriptions <- list(
      "Withheld Columns" = sprintf("%s%s", tagList(
        tags$span(
          style = "font-weight:bold;",
          "The following columns are non-metabolite columns providing meta-information about the data:"
        ),
        tags$ul(lapply(c(
          "sample = Identifies sample name",
          "batch = Identifies batch (large sample sets are separated into batches)",
          "class = Identifies sample type",
          "order = Identifies the order in which samples were injected into the instrument"
        ), function(name) {
          tags$li(name)
        }))
      ), if (isTRUE(p$withhold_cols) && !is.null(p$n_withhold)) {
        tagList(
          tags$span(
            style = "font-weight:bold;",
            "The following columns were withheld from correction:"
          ),
          tags$ul(lapply(d$cleaned$withheld_cols, function(name) {
            tags$li(name)
          }))
        )
      } else {
        ""
      }),
      "Imputation Description" = sprintf("%s%s%s", if (d$imputed$qc_str != "nothing to impute") {
        sprintf("Missing QC values are imputed with %s.<br/>",
                d$imputed$qc_str)
      } else {
        "No missing QC values.<br/>"
      }, if (d$imputed$sam_str != "nothing to impute") {
        sprintf("Missing QC values are imputed with %s.", d$imputed$sam_str)
      } else {
        "No missing sample values."
      }, if (d$imputed$qc_str == "nothing to impute" &&
             d$imputed$sam_str == "nothing to impute") {
        ""
      } else if (p$remove_imputed == TRUE) {
        "<br/>Imputed values are removed after correction."
      }
      else {
        ""
      }),
      "Correction Description" = sprintf(
        "Data was corrected using %s. For each metabolite, this method %s This model regresses peak areas in experimental samples, on an individual metabolite basis, against peak areas in pooled quality control samples.",
        d$corrected$str,
        d$corrected$parameters
      ),
      "Transformation Description" = sprintf("%s <br/> %s", d$transformed$str, if (length(d$transformed$withheld_cols) > 0) {
        tagList(
          tags$span(
            style = "font-weight:bold;",
            "The following columns are withheld from the transformation:"
          ),
          tags$ul(lapply(d$transformed$withheld_cols, function(name) {
            tags$li(name)
          }))
        )
      } else {
        ""
      }),
      "Candidate Outliers" = paste(
        "Possible extreme samples are detected by first grouping samples (QC vs non-QC or by class) and computing RSD.",
        "Metabolites with unstable QC RSD (greater than 30%) are not tested for extreme values.",
        "Robust z\u002Dscores are computed for each value within metabolite by median centering and scaling (MAD, IQR/1.349, SD, or 1) within each group.",
        "Candidate extreme values are non-QC sample\u002Dmetabolite pairs with a z\u002Dscore beyond the threshold of 4 for metabolites with stable QC RSD (<= 20%) or 5 for metabolites with borderline QC RSD (20% < QC RSD <= 30%).",
        "Each candidate is then confirmed with a test chosen by group size: Rosner/ESD for n > 25  (records a strength ratio), otherwise Dixon (if uniquely extreme and 3 <= n <= 30) or Grubbs (if extreme).",
        "Tied or ineligible cases can still be confirmed when the sample's squared Mahalanobis distance is flagged (md_only).",
        "Squared Mahalanois distance is computed in the robust PC score space within each group.",
        "We retain PCs to reach at least 80% variance.",
        "Then a robust covariance is computed using the minimum covariance determinant (MCD), Orthogonalized Gnanadesikan\u002DKettenring (OGK), shrinkage, or classical formula depending on sample size and number of PCs retained.",
        "Note: The confirmed candidates are POSSIBLE extreme values. Futher investigation should be done before removing the metabolite values."
      ),
      "Metabolite Scatter Plots" = sprintf(
        "These plots show metabolites before and after signal drift correction before any transformation is applied. The two metabolites shown above have the largest decrease in sample variation. The change in variation was determined by calculating relative standard deviation (RSD) for each metabolite %s %s%s A full explanation of RSD is in the next section.",
        if (p$rsd_cal == "class_met")
          "grouping by sample class."
        else
          "",
        if (isTRUE(!p$post_cor_filter))
          "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above "
        else
          "",
        if (isTRUE(!p$post_cor_filter))
          sprintf("%s%%.", p$rsd_cutoff)
        else
          ""
      ),
      "RSD Comparison" = sprintf(
        "In these plots, the green indicates RSD decreased after %s, red indicates RSD increased after %s, and gray indicates no change in RSD. For these figures RSD is calculated for each metabolite%s. %s%s <br/> %s ",
        if (p$rsd_compare == "filtered_cor_data")
          "correction"
        else
          "correction and transformation",
        if (p$rsd_compare == "filtered_cor_data")
          "correction"
        else
          "correction and transformation",
        if (p$rsd_cal == "class_met")
          " grouping by sample class"
        else
          "",
        if (isTRUE(!p$post_cor_filter))
          "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above "
        else
          "",
        if (isTRUE(!p$post_cor_filter))
          sprintf("%s%%.", p$rsd_cutoff)
        else
          "",
        if (length(increased_qc) > 0) {
          sprintf(
            "<br/>The following metabolites increased QC RSD after correction: <br/> %s <br/> More investiagtion is needed to determine if these metabolites should be excluded from the data.",
            paste(increased_qc, collapse = ", ")
          )
        } else {
          ""
        }
      ),
      "PCA Comparison" = sprintf(
        "This PCA plot shows both the raw data and %s data colored by %s.",
        if (p$pca_compare == "filtered_cor_data")
          "corrected"
        else
          "corrected and transformed",
        p$color_col
      )
    )
    
    params <- list(
      title = "Metabolomics Data Quality Report",
      notes = p$notes %||% "",
      plots = list(
        "Metabolite Scatter 1" = met1_plot,
        "Metabolite Scatter 2" = met2_plot,
        "RSD Comparison" = rsd_plot,
        "PCA Comparison" = pca_plot,
        "PCA Loading" = pca_loading_plot
      ),
      p = p,
      d = d,
      descriptions = descriptions
    )
    shiny::incProgress(1 / 3, detail = "Saved: report information")
    
    html_out <- rmarkdown::render(
      input         = template,
      output_format = "html_document",
      output_file   = file.path(out_dir, "quality_report.html"),
      params        = params,
      envir         = env,
      quiet         = TRUE
    )
    shiny::incProgress(1 / 3, detail = "Saved: HTML")
    
    list(
      html = normalizePath(html_out, winslash = "/")
    )
  })
}
