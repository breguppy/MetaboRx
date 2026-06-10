#' RSD functions and helpers
#'
#' @keywords internal
#' @noRd

lab_levels <- c("Increased", "No Change", "Decreased")

#' @keywords internal
#' @noRd
color_values <- c(
  "Increased" = "#B22222",
  "No Change" = "gray25",
  "Decreased" = "#234F1E"
)

# ------------------------------------------------------------------------------
# Core utilities
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.get_rsd_data_after <- function(compare_to, p, d) {
  if (compare_to == "filtered_cor_data") {
    title <- "Post-correction Changes"
    sheet_label <- "Corrected"

    if (isTRUE(p$remove_imputed)) {
      df_after <- d$filtered_corrected$df_mv
    } else {
      df_after <- d$filtered_corrected$df_no_mv
    }
  } else {
    title <- "Post-transformation Changes"
    sheet_label <- "Transformed"

    if (isTRUE(p$remove_imputed)) {
      df_after <- d$transformed$df_mv
    } else {
      df_after <- d$transformed$df_no_mv
    }
  }

  list(
    df_before = d$filtered$df,
    df_after = df_after,
    title = title,
    sheet_label = sheet_label
  )
}

#' @keywords internal
#' @noRd
.rsd_change_label <- function(before, after) {
  dplyr::case_when(
    after > before ~ "Increased",
    after < before ~ "Decreased",
    TRUE ~ "No Change"
  )
}

#' @keywords internal
#' @noRd
.rsd_pct_increase <- function(delta) {
  d <- delta[!is.na(delta)]
  if (length(d) == 0L) {
    return(NA_real_)
  }
  100 * mean(d > 0)
}

#' @keywords internal
#' @noRd
.rsd_pct_decrease <- function(delta) {
  d <- delta[!is.na(delta)]
  if (length(d) == 0L) {
    return(NA_real_)
  }
  100 * mean(d < 0)
}

# ------------------------------------------------------------------------------
# Raw RSD calculators
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
metabolite_rsd <- function(df,
                           metadata_cols = c("sample", "batch", "class", "order")) {
  nm <- names(df)
  md_idx <- tolower(nm) %in% tolower(metadata_cols)

  class_col <- nm[tolower(nm) == "class"]
  if (!length(class_col)) {
    stop("Expected a 'class' column (any case).")
  }
  class_col <- class_col[1]

  metab_cols <- nm[!md_idx]
  is_num <- vapply(df[metab_cols], is.numeric, logical(1))
  metab_cols <- metab_cols[is_num]

  if (!length(metab_cols)) {
    stop("No numeric metabolite columns detected.")
  }

  qc_df <- df[df[[class_col]] == "QC", metab_cols, drop = FALSE]
  nonqc_df <- df[df[[class_col]] != "QC", metab_cols, drop = FALSE]

  rsd_fun <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- stats::sd(x, na.rm = TRUE)

    if (!is.finite(mu) || mu == 0) {
      return(NA_real_)
    }

    100 * sigma / mu
  }

  data.frame(
    Metabolite = metab_cols,
    RSD_QC = vapply(qc_df, rsd_fun, numeric(1)),
    RSD_NonQC = vapply(nonqc_df, rsd_fun, numeric(1)),
    check.names = FALSE
  )
}

#' @keywords internal
#' @noRd
class_metabolite_rsd <- function(df,
                                 metadata_cols = c("sample", "batch", "class", "order")) {
  nm <- names(df)
  class_col <- nm[tolower(nm) == "class"]

  if (!length(class_col)) {
    stop("Expected a 'class' column (any case).")
  }
  class_col <- class_col[1]

  meta_idx <- tolower(nm) %in% tolower(metadata_cols)
  metab_cols <- nm[!meta_idx]
  is_num <- vapply(df[metab_cols], is.numeric, logical(1))
  metab_cols <- metab_cols[is_num]

  if (!length(metab_cols)) {
    stop("No numeric metabolite columns detected.")
  }

  long_df <- tidyr::pivot_longer(
    data = df,
    cols = dplyr::all_of(metab_cols),
    names_to = "Metabolite",
    values_to = "Value"
  )

  long_df |>
    dplyr::group_by(.data[[class_col]], Metabolite) |>
    dplyr::summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = stats::sd(Value, na.rm = TRUE),
      RSD = dplyr::if_else(
        is.na(Mean) | Mean == 0,
        NA_real_,
        100 * SD / Mean
      ),
      .groups = "drop"
    ) |>
    dplyr::rename(class = dplyr::all_of(class_col))
}

# ------------------------------------------------------------------------------
# Shared comparison builders
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.compare_metabolite_rsd <- function(rsd_before, rsd_after) {
  df <- dplyr::inner_join(
    dplyr::rename(
      rsd_before,
      QC_before = RSD_QC,
      Samples_before = RSD_NonQC
    ),
    dplyr::rename(
      rsd_after,
      QC_after = RSD_QC,
      Samples_after = RSD_NonQC
    ),
    by = "Metabolite"
  )

  dplyr::bind_rows(
    df |>
      dplyr::transmute(
        analysis = "metabolite",
        Type = "QC",
        unit_id = Metabolite,
        class = NA_character_,
        Metabolite = Metabolite,
        before = QC_before,
        after = QC_after,
        delta = QC_after - QC_before
      ),
    df |>
      dplyr::transmute(
        analysis = "metabolite",
        Type = "Samples",
        unit_id = Metabolite,
        class = NA_character_,
        Metabolite = Metabolite,
        before = Samples_before,
        after = Samples_after,
        delta = Samples_after - Samples_before
      )
  ) |>
    dplyr::mutate(
      change = .rsd_change_label(before, after),
      Type = factor(Type, levels = c("Samples", "QC")),
      change = factor(change, levels = lab_levels)
    )
}

#' @keywords internal
#' @noRd
.compare_class_metabolite_rsd <- function(rsd_before, rsd_after) {
  dplyr::inner_join(
    dplyr::rename(rsd_before, before = RSD),
    dplyr::rename(rsd_after, after = RSD),
    by = c("class", "Metabolite")
  ) |>
    dplyr::mutate(
      analysis = "class_metabolite",
      Type = ifelse(class == "QC", "QC", "Samples"),
      unit_id = paste(class, Metabolite, sep = "::"),
      delta = after - before,
      change = .rsd_change_label(before, after),
      Type = factor(Type, levels = c("Samples", "QC")),
      change = factor(change, levels = lab_levels)
    ) |>
    dplyr::select(
      analysis,
      Type,
      unit_id,
      class,
      Metabolite,
      before,
      after,
      delta,
      change
    )
}

#' @keywords internal
#' @noRd
.summarize_rsd_comparison <- function(compare_df) {
  compare_df |>
    dplyr::group_by(Type) |>
    dplyr::summarise(
      avg_delta = mean(delta, na.rm = TRUE),
      med_delta = stats::median(delta, na.rm = TRUE),
      pct_increase = .rsd_pct_increase(delta),
      pct_decrease = .rsd_pct_decrease(delta),
      .groups = "drop"
    )
}

#' @keywords internal
#' @noRd
.build_rsd_results <- function(df_before, df_after) {
  rsd_before_met <- metabolite_rsd(df_before)
  rsd_after_met <- metabolite_rsd(df_after)

  rsd_before_class <- class_metabolite_rsd(df_before)
  rsd_after_class <- class_metabolite_rsd(df_after)

  compare_met <- .compare_metabolite_rsd(rsd_before_met, rsd_after_met)
  compare_class <- .compare_class_metabolite_rsd(rsd_before_class, rsd_after_class)

  list(
    metabolite = list(
      before = rsd_before_met,
      after = rsd_after_met,
      compare = compare_met,
      stats = .summarize_rsd_comparison(compare_met)
    ),
    class_metabolite = list(
      before = rsd_before_class,
      after = rsd_after_class,
      compare = compare_class,
      stats = .summarize_rsd_comparison(compare_class)
    )
  )
}
