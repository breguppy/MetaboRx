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