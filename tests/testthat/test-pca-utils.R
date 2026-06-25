# test-pca-utils.R

make_pca_pair <- function() {
  before_scores <- data.frame(
    PC1 = c(-1, 0, 1, 2),
    PC2 = c(1, -1, 0, 2),
    sample = paste0("s", 1:4),
    batch = c("b1", "b1", "b2", "b2"),
    class = c("QC", "sample", "QC", "sample"),
    order = 1:4,
    check.names = FALSE
  )

  after_scores <- data.frame(
    PC1 = c(-0.5, 0.2, 0.8, 1.5),
    PC2 = c(0.7, -0.8, 0.1, 1.3),
    sample = paste0("s", 1:4),
    batch = c("b1", "b1", "b2", "b2"),
    class = c("QC", "sample", "QC", "sample"),
    order = 1:4,
    check.names = FALSE
  )

  explained_variance <- data.frame(
    PC = c("PC1", "PC2"),
    explained_variance = c(0.6, 0.3),
    cumulative_explained_variance = c(0.6, 0.9)
  )

  list(
    before = list(scores = before_scores, explained_variance = explained_variance),
    after = list(scores = after_scores, explained_variance = explained_variance)
  )
}

make_large_pca_pair <- function(n = 30L) {
  explained_variance <- data.frame(
    PC = c("PC1", "PC2"),
    explained_variance = c(0.6, 0.3),
    cumulative_explained_variance = c(0.6, 0.9)
  )

  before_scores <- data.frame(
    PC1 = seq(-1, 1, length.out = n),
    PC2 = seq(1, -1, length.out = n),
    sample = paste0("b", seq_len(n)),
    class = rep(c("QC", "sample"), length.out = n),
    many_colors = paste0("c", seq_len(n)),
    many_shapes = paste0("s", seq_len(n)),
    check.names = FALSE
  )

  after_scores <- data.frame(
    PC1 = seq(-0.5, 0.5, length.out = n),
    PC2 = seq(0.5, -0.5, length.out = n),
    sample = paste0("a", seq_len(n)),
    class = rep(c("QC", "sample"), length.out = n),
    many_colors = paste0("c", seq_len(n) + n),
    many_shapes = paste0("s", seq_len(n) + n),
    check.names = FALSE
  )

  list(
    before = list(scores = before_scores, explained_variance = explained_variance),
    after = list(scores = after_scores, explained_variance = explained_variance)
  )
}

test_that("plot_pca_from_result accepts factor-valued shape settings from expand.grid", {
  testthat::skip_if_not_installed("cowplot")

  p <- list(
    color_col = factor("class"),
    shape_col = factor("batch")
  )

  plot <- plot_pca_from_result(
    p = p,
    pca_pair = make_pca_pair(),
    compared_to = "Corrected Data"
  )

  testthat::expect_s3_class(plot, "ggplot")
})

test_that("PCA color choices keep numeric-like columns regardless of group count", {
  meta_df <- data.frame(
    sample = paste0("s", seq_len(30)),
    order_text = as.character(seq_len(30)),
    too_many_groups = paste0("g", seq_len(30)),
    check.names = FALSE
  )

  choices <- .pca_valid_color_choices(meta_df)

  testthat::expect_true("order_text" %in% choices)
  testthat::expect_false("too_many_groups" %in% choices)
})

test_that("PCA shape choices exclude columns with too many groups", {
  meta_df <- data.frame(
    sample = paste0("s", seq_len(30)),
    batch = rep(c("b1", "b2"), length.out = 30),
    too_many_shapes = paste0("s", seq_len(30)),
    check.names = FALSE
  )

  choices <- .pca_valid_shape_choices(meta_df)

  testthat::expect_true("batch" %in% choices)
  testthat::expect_false("too_many_shapes" %in% choices)
})

test_that("plot_pca_from_result uses gradient for numeric-like color columns", {
  testthat::skip_if_not_installed("cowplot")

  pca_pair <- make_pca_pair()
  pca_pair$before$scores$order_text <- as.character(pca_pair$before$scores$order)
  pca_pair$after$scores$order_text <- as.character(pca_pair$after$scores$order)

  plot <- plot_pca_from_result(
    p = list(color_col = "order_text", shape_col = "none"),
    pca_pair = pca_pair,
    compared_to = "Corrected Data"
  )

  testthat::expect_s3_class(plot, "ggplot")
})

test_that("plot_pca_from_result errors for forced over-limit color and shape groups", {
  pca_pair <- make_large_pca_pair()

  testthat::expect_error(
    plot_pca_from_result(
      p = list(color_col = "many_colors", shape_col = "none"),
      pca_pair = pca_pair,
      compared_to = "Corrected Data"
    ),
    "Too many groups for PCA color palette"
  )

  testthat::expect_error(
    plot_pca_from_result(
      p = list(color_col = "class", shape_col = "many_shapes"),
      pca_pair = pca_pair,
      compared_to = "Corrected Data"
    ),
    "Too many groups for point shapes"
  )
})

test_that("plot_pca_from_result accepts explicit no-shape option", {
  testthat::skip_if_not_installed("cowplot")

  plot <- plot_pca_from_result(
    p = list(color_col = "class", shape_col = "none"),
    pca_pair = make_pca_pair(),
    compared_to = "Corrected Data"
  )

  testthat::expect_s3_class(plot, "ggplot")
})

test_that("PCA UI offers no-shape choice and filters over-limit columns", {
  meta_df <- data.frame(
    sample = paste0("s", seq_len(30)),
    batch = rep(c("b1", "b2"), length.out = 30),
    order_text = as.character(seq_len(30)),
    too_many_groups = paste0("g", seq_len(30)),
    too_many_shapes = paste0("s", seq_len(30)),
    check.names = FALSE
  )

  rendered <- htmltools::renderTags(ui_pca_eval(meta_df, ns = identity))$html

  testthat::expect_match(rendered, "value=\"none\"", fixed = TRUE)
  testthat::expect_match(rendered, "No shapes", fixed = TRUE)
  testthat::expect_match(rendered, "value=\"order_text\"", fixed = TRUE)
  testthat::expect_false(grepl("value=\"too_many_groups\"", rendered, fixed = TRUE))
  testthat::expect_false(grepl("value=\"too_many_shapes\"", rendered, fixed = TRUE))
})

make_pca_export_data <- function(drop_corrected_met_c = FALSE) {
  raw_df <- data.frame(
    sample = paste0("s", seq_len(6)),
    batch = rep(c("b1", "b2"), each = 3),
    class = c("QC", "sample", "sample", "QC", "sample", "sample"),
    order = seq_len(6),
    met_a = c(10, 11, 12, 10.5, 11.5, 12.5),
    met_b = c(20, 19, 21, 20.5, 19.5, 21.5),
    met_c = c(30, 31, 29, 30.5, 31.5, 29.5),
    check.names = FALSE
  )

  corrected_df <- dplyr::mutate(
    raw_df,
    met_a = met_a / mean(met_a),
    met_b = met_b / mean(met_b),
    met_c = met_c / mean(met_c)
  )

  if (isTRUE(drop_corrected_met_c)) {
    corrected_df <- corrected_df[, setdiff(names(corrected_df), "met_c"), drop = FALSE]
  }

  list(
    filtered = list(df = raw_df),
    filtered_corrected = list(df_no_mv = corrected_df, df_mv = corrected_df),
    transformed = list(df_no_mv = corrected_df, df_mv = corrected_df)
  )
}

test_that("PCA export mapping reuses paired results when metabolite columns match", {
  p <- list(remove_imputed = FALSE, transform = "none")
  d <- make_pca_export_data()
  pair <- compute_pca_pair(
    before = d$filtered$df,
    after = d$filtered_corrected$df_no_mv,
    p = p
  )

  reusable <- .pca_export_results_from_pairs(
    p = p,
    datasets = get_pca_export_datasets(p, d),
    pca_pairs = list(filtered_cor_data = pair)
  )

  testthat::expect_identical(reusable$raw_data, pair$before)
  testthat::expect_identical(reusable$corrected_data, pair$after)
})

test_that("PCA export mapping falls back when reuse would change metabolite columns", {
  p <- list(remove_imputed = FALSE, transform = "none")
  d <- make_pca_export_data(drop_corrected_met_c = TRUE)
  pair <- compute_pca_pair(
    before = d$filtered$df,
    after = d$filtered_corrected$df_no_mv,
    p = p
  )

  with_reuse <- compute_all_pca_export_results(
    p = p,
    d = d,
    pca_pairs = list(filtered_cor_data = pair)
  )
  without_reuse <- compute_all_pca_export_results(p = p, d = d)

  testthat::expect_identical(with_reuse$raw_data$metab_cols, c("met_a", "met_b", "met_c"))
  testthat::expect_identical(with_reuse$corrected_data$metab_cols, c("met_a", "met_b"))
  testthat::expect_equal(
    with_reuse$raw_data$explained_variance,
    without_reuse$raw_data$explained_variance,
    tolerance = 1e-12
  )
})
