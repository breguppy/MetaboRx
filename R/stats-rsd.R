#' RSD functions
#'
#' @keywords internal
#' @noRd
metabolite_rsd <- function(df,
                           metadata_cols = c("sample", "batch", "class", "order")) {
  nm <- names(df)
  md_idx <- tolower(nm) %in% tolower(metadata_cols)
  class_col <- nm[tolower(nm) == "class"]
  if (!length(class_col))
    stop("Expected a 'class' column (any case).")
  class_col <- class_col[1]
  metab_cols <- nm[!md_idx]
  is_num <- vapply(df[metab_cols], is.numeric, logical(1))
  metab_cols <- metab_cols[is_num]
  if (!length(metab_cols))
    stop("No numeric metabolite columns detected.")
  
  # Separate QC and non-QC samples
  qc_df <- df[df[[class_col]] == "QC", metab_cols, drop = FALSE]
  nonqc_df <- df[df[[class_col]] != "QC", metab_cols, drop = FALSE]
  
  # RSD function ignoring NA
  rsd_fun <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- stats::sd(x, na.rm = TRUE)
    if (mu == 0 || !is.finite(mu))
      return(NA_real_)
    return(100 * sigma / mu)
  }
  
  data.frame(
    Metabolite = metab_cols,
    RSD_QC     = vapply(qc_df, rsd_fun, numeric(1)),
    RSD_NonQC  = vapply(nonqc_df, rsd_fun, numeric(1)),
    check.names = FALSE
  )
}

#' @keywords internal
#' @noRd
class_metabolite_rsd <- function(df,
                                 metadata_cols = c("sample", "batch", "class", "order")) {
  # Get metabolite columns
  metab_cols <- setdiff(names(df), metadata_cols)
  
  long_df <- df %>%
    pivot_longer(
      cols = all_of(metab_cols),
      names_to = "Metabolite",
      values_to = "Value"
    )
  
  # Group by class and metabolite and compute RSD
  rsd_df <- long_df %>%
    group_by(class, Metabolite) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      RSD = ifelse(Mean == 0, NA, SD / Mean * 100),
      .groups = "drop"
    )
  
  return(rsd_df)
}

#' @keywords internal
#' @noRd
#' @keywords internal
#' @noRd
delta_rsd_stats <- function(rsdBefore, rsdAfter) {
  mean_na   <- function(x) mean(x, na.rm = TRUE)
  median_na <- function(x) stats::median(x, na.rm = TRUE)
  
  pct_increase <- function(delta) {
    d <- delta[!is.na(delta)]
    if (length(d) == 0L) return(NA_real_)
    100 * mean(d > 0)
  }
  
  pct_decrease <- function(delta) {
    d <- delta[!is.na(delta)]
    if (length(d) == 0L) return(NA_real_)
    100 * mean(d < 0)
  }
  
  if (all(c("Metabolite", "RSD_QC", "RSD_NonQC") %in% names(rsdBefore)) &&
      all(c("Metabolite", "RSD_QC", "RSD_NonQC") %in% names(rsdAfter))) {
    
    df <- dplyr::inner_join(
      dplyr::rename(rsdBefore, RSD_QC_b = RSD_QC, RSD_NonQC_b = RSD_NonQC),
      dplyr::rename(rsdAfter,  RSD_QC_a = RSD_QC, RSD_NonQC_a = RSD_NonQC),
      by = "Metabolite"
    )
    
    delta_qc <- df$RSD_QC_a    - df$RSD_QC_b
    delta_s  <- df$RSD_NonQC_a - df$RSD_NonQC_b
    
  } else if (all(c("class", "Metabolite", "RSD") %in% names(rsdBefore)) &&
             all(c("class", "Metabolite", "RSD") %in% names(rsdAfter))) {
    
    df <- dplyr::inner_join(
      dplyr::rename(rsdBefore, RSD_b = RSD),
      dplyr::rename(rsdAfter,  RSD_a = RSD),
      by = c("class", "Metabolite")
    ) |>
      dplyr::mutate(delta = RSD_a - RSD_b)
    
    delta_qc <- df |>
      dplyr::filter(class == "QC") |>
      dplyr::pull(delta)
    
    # For samples: average delta across non-QC classes within metabolite (as you already do)
    delta_s <- df |>
      dplyr::filter(class != "QC") |>
      dplyr::group_by(Metabolite) |>
      dplyr::summarise(delta_s = mean(delta, na.rm = TRUE), .groups = "drop") |>
      dplyr::pull(delta_s)
    
  } else {
    stop("Unrecognized input schema for rsdBefore/rsdAfter.")
  }
  
  list(
    avg_delta_qc          = mean_na(delta_qc),
    med_delta_qc          = median_na(delta_qc),
    avg_delta_sample      = mean_na(delta_s),
    med_delta_sample      = median_na(delta_s),
    pct_increase_qc       = pct_increase(delta_qc),
    pct_decrease_qc       = pct_decrease(delta_qc),
    pct_increase_sample   = pct_increase(delta_s),
    pct_decrease_sample   = pct_decrease(delta_s)
  )
}
