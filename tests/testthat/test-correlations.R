test_that("compute_pairwise_metabolite_correlations returns sorted finite pairs", {
  df <- data.frame(
    sample = paste0("s", 1:5),
    batch = 1L,
    class = c("QC", "sample", "sample", "QC", "sample"),
    order = 1:5,
    A = c(1, 2, 3, 4, 5),
    B = c(2, 4, 6, 8, 10),
    C = c(5, 4, 3, 2, 1),
    D = c(1, NA, 1, NA, 1),
    stringsAsFactors = FALSE
  )

  out <- compute_pairwise_metabolite_correlations(
    df = df,
    cols = c("A", "B", "C", "D"),
    min_complete = 3L
  )

  expect_named(out, c("col1", "col2", "cor", "n_complete"))
  expect_equal(out$col1[[1]], "A")
  expect_equal(out$col2[[1]], "B")
  expect_equal(out$cor[[1]], 1)
  expect_true(all(abs(out$cor) == sort(abs(out$cor), decreasing = TRUE)))
  expect_true(all(out$n_complete >= 3L))
})

test_that("filter_correlation_pairs_by_range preserves descending absolute order", {
  corr_df <- data.frame(
    col1 = c("A", "A", "B", "C"),
    col2 = c("B", "C", "C", "D"),
    cor = c(0.95, -0.99, 0.90, 0.97),
    n_complete = c(5L, 5L, 5L, 5L),
    stringsAsFactors = FALSE
  )

  out <- filter_correlation_pairs_by_range(corr_df, c(0.94, 1))

  expect_equal(out$cor, c(0.97, 0.95))
  expect_equal(out$col1, c("C", "A"))
})
