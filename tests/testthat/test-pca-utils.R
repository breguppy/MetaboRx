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
