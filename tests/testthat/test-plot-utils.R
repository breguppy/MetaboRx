test_that("RSD facet swatches remain portable to the base PDF device", {
  swatch <- circle("#336699")

  expect_match(swatch, "color:#336699", fixed = TRUE)
  expect_false(grepl("[^ -~]", swatch))
  expect_false(grepl("&#[0-9]+;", swatch))
})
