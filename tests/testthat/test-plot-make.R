make_rsd_export_df <- function() {
  data.frame(
    sample = paste0("S", seq_len(8)),
    batch = rep(c("B1", "B2"), each = 4),
    class = c("QC", "QC", "Case", "Control", "QC", "QC", "Case", "Control"),
    order = seq_len(8),
    met_a = c(10, 12, 30, 35, 11, 13, 31, 37),
    met_b = c(50, 48, 70, 72, 49, 47, 69, 73),
    check.names = FALSE
  )
}

make_rsd_export_data <- function() {
  df_before <- make_rsd_export_df()
  df_corrected <- dplyr::mutate(
    df_before,
    met_a = met_a * c(1, 0.96, 0.75, 0.72, 1.02, 0.98, 0.76, 0.74),
    met_b = met_b * c(1, 1.01, 0.85, 0.82, 0.99, 1.02, 0.86, 0.84)
  )
  df_transformed <- dplyr::mutate(
    df_corrected,
    met_a = log1p(met_a),
    met_b = log1p(met_b)
  )

  list(
    filtered = list(df = df_before),
    filtered_corrected = list(df_no_mv = df_corrected, df_mv = df_corrected),
    transformed = list(df_no_mv = df_transformed, df_mv = df_transformed)
  )
}

make_rsd_export_params <- function(transform = "none") {
  list(
    remove_imputed = FALSE,
    transform = transform
  )
}

with_stubbed_rsd_scatter_plot <- function(code) {
  ns <- asNamespace("QCcorrection")
  old_mk_plot <- get("mk_plot", envir = ns)
  was_locked <- bindingIsLocked("mk_plot", ns)

  if (was_locked) {
    unlockBinding("mk_plot", ns)
  }
  assign(
    "mk_plot",
    function(d_all, x, y, facet_labs, compared_to) {
      ggplot2::ggplot(d_all, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
        ggplot2::labs(
          title = paste("Comparison of RSD Before and After", compared_to),
          x = "RSD (%) Before",
          y = "RSD (%) After"
        )
    },
    envir = ns
  )
  if (was_locked) {
    lockBinding("mk_plot", ns)
  }

  on.exit(
    {
      if (bindingIsLocked("mk_plot", ns)) {
        unlockBinding("mk_plot", ns)
      }
      assign("mk_plot", old_mk_plot, envir = ns)
      if (was_locked) {
        lockBinding("mk_plot", ns)
      }
    },
    add = TRUE
  )

  force(code)
}

with_counted_rsd_results <- function(counter, code) {
  ns <- asNamespace("QCcorrection")
  old_build <- get(".build_rsd_results", envir = ns)
  was_locked <- bindingIsLocked(".build_rsd_results", ns)

  if (was_locked) {
    unlockBinding(".build_rsd_results", ns)
  }
  assign(
    ".build_rsd_results",
    function(df_before, df_after) {
      counter(counter() + 1L)
      old_build(df_before, df_after)
    },
    envir = ns
  )
  if (was_locked) {
    lockBinding(".build_rsd_results", ns)
  }

  on.exit(
    {
      if (bindingIsLocked(".build_rsd_results", ns)) {
        unlockBinding(".build_rsd_results", ns)
      }
      assign(".build_rsd_results", old_build, envir = ns)
      if (was_locked) {
        lockBinding(".build_rsd_results", ns)
      }
    },
    add = TRUE
  )

  force(code)
}

test_that("make_all_rsd_plots creates distinct correction-only export specs", {
  res <- with_stubbed_rsd_scatter_plot(
    make_all_rsd_plots(
      make_rsd_export_params(transform = "none"),
      make_rsd_export_data()
    )
  )

  testthat::expect_length(res[["rsd_plots"]], 4L)
  testthat::expect_setequal(
    res[["plot_names"]],
    c(
      "rsd_dist_filtered_cor_data_class_met",
      "rsd_scatter_filtered_cor_data_class_met",
      "rsd_dist_filtered_cor_data_met",
      "rsd_scatter_filtered_cor_data_met"
    )
  )
  testthat::expect_false(any(grepl("class-met", res[["plot_names"]], fixed = TRUE)))

  y_labels <- purrr::map_chr(res[["rsd_plots"]], function(plot) plot[["labels"]][["y"]])
  names(y_labels) <- res[["plot_names"]]

  testthat::expect_true(all(y_labels[grepl("rsd_dist", names(y_labels))] == "Density"))
  testthat::expect_true(all(y_labels[grepl("rsd_scatter", names(y_labels))] == "RSD (%) After"))
})

test_that("make_all_rsd_plots includes transformed RSD export specs when enabled", {
  res <- with_stubbed_rsd_scatter_plot(
    make_all_rsd_plots(
      make_rsd_export_params(transform = "log"),
      make_rsd_export_data()
    )
  )

  testthat::expect_length(res[["rsd_plots"]], 8L)
  testthat::expect_true(any(grepl("filtered_cor_data", res[["plot_names"]], fixed = TRUE)))
  testthat::expect_true(any(grepl("transformed_cor_data", res[["plot_names"]], fixed = TRUE)))
})

test_that("make_all_rsd_plots reuses RSD results for each comparison target", {
  counter <- local({
    count <- 0L
    function(value) {
      if (!missing(value)) {
        count <<- value
      }
      count
    }
  })

  res <- with_stubbed_rsd_scatter_plot(
    with_counted_rsd_results(
      counter,
      make_all_rsd_plots(
        make_rsd_export_params(transform = "none"),
        make_rsd_export_data()
      )
    )
  )

  testthat::expect_length(res[["rsd_plots"]], 4L)
  testthat::expect_equal(counter(), 1L)

  counter(0L)
  res <- with_stubbed_rsd_scatter_plot(
    with_counted_rsd_results(
      counter,
      make_all_rsd_plots(
        make_rsd_export_params(transform = "log"),
        make_rsd_export_data()
      )
    )
  )

  testthat::expect_length(res[["rsd_plots"]], 8L)
  testthat::expect_equal(counter(), 2L)
})

test_that("RSD plots can be built from precomputed RSD results", {
  d <- make_rsd_export_data()
  rsd_results <- .build_rsd_results(d$filtered$df, d$filtered_corrected$df_no_mv)

  scatter_plot <- plot_rsd_comparison_from_results(rsd_results, "Correction")
  dist_plot <- plot_met_rsd_distributions_from_results(rsd_results, "Correction")

  testthat::expect_s3_class(scatter_plot, "ggplot")
  testthat::expect_s3_class(dist_plot, "ggplot")
  testthat::expect_equal(scatter_plot[["labels"]][["x"]], "RSD (%) Before")
  testthat::expect_equal(scatter_plot[["labels"]][["y"]], "RSD (%) After")
  testthat::expect_equal(dist_plot[["labels"]][["y"]], "Density")
})

test_that("make_rsd_plot validates RSD plot options", {
  d <- make_rsd_export_data()
  base_params <- make_rsd_export_params(transform = "none")

  legacy_class_params <- utils::modifyList(
    base_params,
    list(
      rsd_compare = "filtered_cor_data",
      rsd_cal = "class-met",
      rsd_plot_type = "dist"
    )
  )
  testthat::expect_s3_class(make_rsd_plot(legacy_class_params, d), "ggplot")

  bad_plot_type <- utils::modifyList(legacy_class_params, list(rsd_plot_type = "box"))
  testthat::expect_error(make_rsd_plot(bad_plot_type, d), "rsd_plot_type")

  bad_rsd_cal <- utils::modifyList(legacy_class_params, list(rsd_cal = "by_batch"))
  testthat::expect_error(make_rsd_plot(bad_rsd_cal, d), "rsd_cal")

  bad_compare <- utils::modifyList(legacy_class_params, list(rsd_compare = "QC vs Sample"))
  testthat::expect_error(make_rsd_plot(bad_compare, d), "rsd_compare")
})
