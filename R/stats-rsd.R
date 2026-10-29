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
#' Calculate RSD changes with numerical tolerance
#'
#' Differences within tolerance are treated as zero.
#' Nonfinite before/after pairs are excluded using NA.
#'
#' @keywords internal
#' @noRd
.rsd_change_delta <- function(
    before,
    after,
    atol = 1e-8,
    rtol = 1e-10
) {
  if (!is.numeric(before) || !is.numeric(after)) {
    stop("before and after must be numeric vectors.")
  }
  
  if (length(before) != length(after)) {
    stop("before and after must have the same length.")
  }
  
  if (!is.numeric(atol) ||
      length(atol) != 1L ||
      !is.finite(atol) ||
      atol < 0) {
    stop("atol must be a finite, non-negative numeric scalar.")
  }
  
  if (!is.numeric(rtol) ||
      length(rtol) != 1L ||
      !is.finite(rtol) ||
      rtol < 0) {
    stop("rtol must be a finite, non-negative numeric scalar.")
  }
  
  delta <- rep(NA_real_, length(before))
  valid <- is.finite(before) & is.finite(after)
  
  if (any(valid)) {
    differences <- after[valid] - before[valid]
    tolerance <- atol +
      rtol * pmax(abs(before[valid]), abs(after[valid]))
    
    # Exclude differences that overflow despite finite inputs.
    finite_difference <- is.finite(differences)
    unchanged <- finite_difference &
      abs(differences) <= tolerance
    
    differences[unchanged] <- 0
    differences[!finite_difference] <- NA_real_
    
    delta[valid] <- differences
  }
  
  delta
}

#' @keywords internal
#' @noRd
.rsd_change_label <- function(
    before,
    after,
    atol = 1e-8,
    rtol = 1e-10
) {
  delta <- .rsd_change_delta(
    before = before,
    after = after,
    atol = atol,
    rtol = rtol
  )
  
  dplyr::case_when(
    is.na(delta) ~ NA_character_,
    delta > 0 ~ "Increased",
    delta < 0 ~ "Decreased",
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

  rows_by_class <- split(seq_len(nrow(df)), df[[class_col]])
  rsd_by_class <- vector("list", length(rows_by_class))
  names(rsd_by_class) <- names(rows_by_class)

  for (i in seq_along(rows_by_class)) {
    idx <- rows_by_class[[i]]
    metab_df <- df[idx, metab_cols, drop = FALSE]

    means <- vapply(metab_df, mean, numeric(1L), na.rm = TRUE)
    sds <- vapply(metab_df, stats::sd, numeric(1L), na.rm = TRUE)
    rsd <- dplyr::if_else(
      is.na(means) | means == 0,
      NA_real_,
      100 * sds / means
    )

    rsd_by_class[[i]] <- data.frame(
      class = df[[class_col]][idx[[1L]]],
      Metabolite = metab_cols,
      Mean = unname(means),
      SD = unname(sds),
      RSD = unname(rsd),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rsd_by_class)
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
        after = QC_after
      ),
    df |>
      dplyr::transmute(
        analysis = "metabolite",
        Type = "Samples",
        unit_id = Metabolite,
        class = NA_character_,
        Metabolite = Metabolite,
        before = Samples_before,
        after = Samples_after
      )
  ) |>
    dplyr::mutate(
      delta = .rsd_change_delta(before, after),
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
      delta = .rsd_change_delta(before, after),
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
    dplyr::mutate(
      delta = .rsd_change_delta(before, after)
    ) |>
    dplyr::group_by(Type) |>
    dplyr::summarise(
      avg_delta = if (all(is.na(delta))) {
        NA_real_
      } else {
        mean(delta, na.rm = TRUE)
      },
      med_delta = if (all(is.na(delta))) {
        NA_real_
      } else {
        stats::median(delta, na.rm = TRUE)
      },
      pct_increase = .rsd_pct_increase(delta),
      pct_decrease = .rsd_pct_decrease(delta),
      pct_unchanged = if (all(is.na(delta))) {
        NA_real_
      } else {
        100 * mean(delta == 0, na.rm = TRUE)
      },
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

#' @keywords internal
#' @noRd
.summarize_rsd_changes <- function(
    rsd_before,
    rsd_after,
    atol = 1e-8,
    rtol = 1e-10
) {
  # This call validates paired inputs and tolerances, even for empty inputs.
  delta <- .rsd_change_delta(
    before = rsd_before,
    after = rsd_after,
    atol = atol,
    rtol = rtol
  )
  
  delta <- delta[is.finite(delta)]
  
  if (!length(delta)) {
    return(data.frame(
      med_delta = NA_real_,
      pct_increase = NA_real_,
      pct_decrease = NA_real_,
      pct_unchanged = NA_real_
    ))
  }
  
  data.frame(
    med_delta = stats::median(delta),
    pct_increase = 100 * mean(delta > 0),
    pct_decrease = 100 * mean(delta < 0),
    pct_unchanged = 100 * mean(delta == 0)
  )
}