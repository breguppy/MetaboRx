#' Group stats functions
#'
#' @keywords internal
#' @noRd
group_stats <- function(df) {
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  grouped_df <- df |>
    dplyr::rename(Group = dplyr::all_of("class")) |>
    dplyr::select(-dplyr::any_of(c("batch", "order")))

  summarize_group <- function(group_df, label, fn) {
    group_df |>
      dplyr::summarise(dplyr::across(dplyr::all_of(metab_cols), fn)) |>
      dplyr::mutate(` ` = label)
  }

  group_dfs <- list()
  group_stats_dfs <- list()
  group_names <- unique(grouped_df[["Group"]])

  for (group_name in group_names) {
    group_df <- grouped_df |>
      dplyr::filter(.data$Group == group_name)
    group_dfs[[group_name]] <- group_df

    group_stats_df <- dplyr::bind_rows(
      summarize_group(group_df, "Mean", ~ mean(.x, na.rm = TRUE)),
      summarize_group(group_df, "SE", ~ stats::sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
      summarize_group(group_df, "CV", ~ stats::sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE))
    ) |>
      dplyr::select(` `, dplyr::all_of(metab_cols))
    group_stats_dfs[[group_name]] <- group_stats_df
  }

  return(list(group_dfs = group_dfs, group_stats_dfs = group_stats_dfs))
}

fold_changes <- function(df, control_mean) {
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))

  if (nrow(control_mean) != 1) {
    stop("control_mean must contain exactly one row with metabolite means.")
  }

  control_values <- unlist(control_mean[, metab_cols, drop = FALSE], use.names = TRUE)

  df |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(metab_cols),
        ~ .x / control_values[dplyr::cur_column()]
      )
    )
}
