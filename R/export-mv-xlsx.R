#' exports missing value report as excel file
#'
#' @keywords internal
#' @noRd
export_mv_xlsx <- function(p, d, file = NULL) {
  cleaned_df <- d$cleaned$df

  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  styles <- .xlsx_export_styles()

  shiny::withProgress(message = "Creating missing value summary...", value = 0, {
    meta_cols <- c("sample", "batch", "class", "order")
    metab_cols <- setdiff(names(cleaned_df), meta_cols)
    n_metabs <- length(metab_cols)

    missing_mat <- if (n_metabs == 0L) {
      matrix(FALSE, nrow = nrow(cleaned_df), ncol = 0L)
    } else {
      vapply(
        cleaned_df[, metab_cols, drop = FALSE],
        function(col) is.na(col) | col <= 0,
        logical(nrow(cleaned_df))
      )
    }
    if (is.null(dim(missing_mat))) {
      missing_mat <- matrix(missing_mat, nrow = nrow(cleaned_df))
    }
    colnames(missing_mat) <- metab_cols

    count_missing_by_metabolite <- function(row_idx) {
      n_samples <- length(row_idx)
      missing_count <- if (n_metabs == 0L) {
        integer(0)
      } else {
        colSums(missing_mat[row_idx, , drop = FALSE])
      }

      tibble::tibble(
        metabolite = metab_cols,
        missing_count = as.integer(missing_count),
        missing_pct = if (n_samples == 0L) {
          NA_real_
        } else {
          (missing_count / n_samples) * 100
        }
      ) |>
        dplyr::arrange(dplyr::desc(.data$missing_pct))
    }

    sample_idx <- which(cleaned_df[["class"]] != "QC")
    qc_idx <- which(cleaned_df[["class"]] == "QC")

    sample_renamed <- count_missing_by_metabolite(sample_idx) |>
      dplyr::transmute(
        metabolite = .data$metabolite,
        sample_missing_count = .data$missing_count,
        sample_missing_pct = .data$missing_pct
      )

    qc_renamed <- count_missing_by_metabolite(qc_idx) |>
      dplyr::transmute(
        metabolite = .data$metabolite,
        qc_missing_count = .data$missing_count,
        qc_missing_pct = .data$missing_pct
      )

    metab_mv <- sample_renamed |>
      dplyr::inner_join(qc_renamed, by = "metabolite") |>
      dplyr::filter(!(.data$sample_missing_pct == 0 & .data$qc_missing_pct == 0)) |>
      dplyr::arrange(dplyr::desc(.data$sample_missing_pct))

    .xlsx_write_described_sheet(
      wb,
      sheet_name = "Metabolite",
      description = paste(
        "Tab Metabolite. Missing value counts (missing_count) and percentages (missing_pct)",
        "per metabolite are listed here for samples and QC samples. Missing values are",
        "defined as NA or <= 0. If a metabolite is not listed here, it did not have any missing values."
      ),
      x = metab_mv,
      merge_cols = 6,
      styles = styles
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by metabolite")

    row_missing_count <- rowSums(missing_mat)
    sample_mv <- tibble::tibble(
      sample = cleaned_df[["sample"]],
      missing_count = row_missing_count,
      missing_pct = if (n_metabs == 0L) {
        NA_real_
      } else {
        (row_missing_count / n_metabs) * 100
      }
    ) |>
      dplyr::filter(.data$missing_pct > 0) |>
      dplyr::arrange(dplyr::desc(.data$missing_pct))

    .xlsx_write_described_sheet(
      wb,
      sheet_name = "Sample",
      description = paste(
        "Tab Sample. Missing value counts (missing_count) and percentages (missing_pct)",
        "per sample are listed here. Missing values are defined as NA or <= 0.",
        "If a sample is not listed here, it did not have any missing values."
      ),
      x = sample_mv,
      merge_cols = 6,
      styles = styles
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by sample")

    summarize_missing_by <- function(group_col) {
      tibble::tibble(
        group_value = cleaned_df[[group_col]],
        row_missing_count = row_missing_count
      ) |>
        dplyr::group_by(.data$group_value) |>
        dplyr::summarise(
          n_samples = dplyr::n(),
          missing_count = sum(.data$row_missing_count),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          total_values = .data$n_samples * n_metabs,
          missing_pct = dplyr::if_else(
            .data$total_values == 0,
            NA_real_,
            (.data$missing_count / .data$total_values) * 100
          )
        ) |>
        dplyr::transmute(
          "{group_col}" := .data$group_value,
          missing_count = .data$missing_count,
          missing_pct = .data$missing_pct
        ) |>
        dplyr::filter(.data$missing_pct > 0) |>
        dplyr::arrange(dplyr::desc(.data$missing_pct))
    }

    class_mv <- summarize_missing_by("class")

    .xlsx_write_described_sheet(
      wb,
      sheet_name = "Class",
      description = paste(
        "Tab Class. Missing value counts (missing_count) and percentages (missing_pct)",
        "per sample class are listed here. Missing values are defined as NA or <= 0.",
        "If a sample class is not listed here, it did not have any missing values."
      ),
      x = class_mv,
      merge_cols = 6,
      styles = styles
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by class")

    batch_mv <- summarize_missing_by("batch")

    .xlsx_write_described_sheet(
      wb,
      sheet_name = "Batch",
      description = paste(
        "Tab Batch. Missing value counts (missing_count) and percentages (missing_pct)",
        "per batch are listed here. Missing values are defined as NA or <= 0.",
        "If a batch is not listed here, it did not have any missing values."
      ),
      x = batch_mv,
      merge_cols = 6,
      styles = styles
    )
    shiny::incProgress(1 / 5, detail = "Saved: missing values by batch")

    class_metab_missing_summary <- if (n_metabs == 0L) {
      data.frame(
        class = character(0),
        metabolite = character(0),
        n_rows_in_class = integer(0),
        missing_count = integer(0),
        missing_pct = numeric(0),
        all_missing = logical(0)
      )
    } else {
      class_summaries <- lapply(sort(unique(cleaned_df[["class"]])), function(class_value) {
        idx <- which(cleaned_df[["class"]] == class_value)
        n_in_class <- length(idx)
        missing_count <- colSums(missing_mat[idx, , drop = FALSE])

        data.frame(
          class = class_value,
          metabolite = metab_cols,
          n_rows_in_class = n_in_class,
          missing_count = as.integer(missing_count),
          missing_pct = (missing_count / n_in_class) * 100,
          all_missing = missing_count == n_in_class,
          row.names = NULL,
          check.names = FALSE
        )
      })

      dplyr::bind_rows(class_summaries)
    }

    if (nrow(class_metab_missing_summary) > 0L) {
      class_metab_missing_export <- class_metab_missing_summary |>
        dplyr::arrange(
          dplyr::desc(.data$all_missing),
          dplyr::desc(.data$missing_pct),
          .data$class,
          .data$metabolite
        )

      .xlsx_write_described_sheet(
        wb,
        sheet_name = "Class-Met Missing",
        description = paste(
          "Tab Class-Met Missing. Missing value counts and percentages are listed",
          "for each class-metabolite pair. Missing values are defined as NA or <= 0.",
          "The column all_missing indicates whether all values for that metabolite",
          "within the class are missing."
        ),
        x = class_metab_missing_export,
        merge_cols = 6,
        styles = styles,
        note_row_height = 75
      )
    }

    shiny::incProgress(1 / 5, detail = "Saved: missing values by class and metabolite")
  })

  .xlsx_save_or_return(wb, file)
}
