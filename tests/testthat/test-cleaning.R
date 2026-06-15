test_that("clean_data basic cleaning and outputs", {
  # toy input (unsorted; QC should start and end after sort + normalization)
  df <- data.frame(
    SampleID  = paste0("s", 1:6),
    BatchID   = c(1, 1, 1, 1, 1, 1),
    Type      = c(NA, "qc", "qc", "Qc", "sample", "QC"),
    Injection = c(3, 1, 6, 2, 5, 4),
    met1      = c("1", "a", "0", "5", "3", NA),
    met2      = c(2, "7", "foo", 0, NA, 1),
    met3      = c("0", "0", "0", "0", "0", "0"),
    note      = c("x", "y", "z", "w", "v", "u"),
    stringsAsFactors = FALSE
  )
  
  out <- clean_data(
    df,
    sample = "SampleID",
    batch  = "BatchID",
    class  = "Type",
    order  = "Injection",
    withheld_cols = c("note")
  )
  
  # structure
  expect_type(out, "list")
  expect_named(out, c(
    "df",
    "meta_df",
    "replacement_counts",
    "withheld_cols",
    "non_numeric_cols",
    "all_missing_zero_qc_cols",
    "duplicate_mets",
    "duplicate_col_names",
    "blank_df"
  ))
  expect_equal(out$withheld_cols, "note")
  expect_equal(out$all_missing_zero_qc_cols, "met3")
  
  # columns renamed, withheld removed, order applied
  expect_setequal(names(out$df),
                  c("sample", "batch", "class", "order", "met1", "met2"))
  expect_true(is.unsorted(df$Injection))
  expect_equal(out$df$order, sort(out$df$order))
  
  # class normalization to "QC"
  expect_true(all(out$df$class %in% c("QC", "sample")))
  expect_identical(out$df$class[1], "QC")
  expect_identical(out$df$class[nrow(out$df)], "QC")
  
  # numeric coercion + zero→NA
  expect_true(all(vapply(out$df[c("met1", "met2")], is.numeric, TRUE)))
  expect_true(any(is.na(out$df$met1)))
  expect_true(any(is.na(out$df$met2)))
  
  # replacement counts
  rc <- out$replacement_counts
  # met1: one non-numeric ("a"), one zero ("0")
  expect_equal(rc$non_numeric_replaced[rc$metabolite == "met1"], 1)
  expect_equal(unname(rc$zero_replaced[rc$metabolite == "met1"]), 1)
  # met2: one non-numeric ("foo"), one zero
  expect_equal(rc$non_numeric_replaced[rc$metabolite == "met2"], 1)
  expect_equal(unname(rc$zero_replaced[rc$metabolite == "met2"]), 1)
})

test_that("clean_data errors when first sample after sort is not QC", {
  df <- data.frame(
    SampleID  = c("s1", "s2", "s3"),
    BatchID   = 1,
    Type      = c("sample", "QC", "QC"),
    Injection = c(1, 2, 3),
    met1      = c(1, 2, 3)
  )
  expect_error(
    clean_data(
      df,
      "SampleID",
      "BatchID",
      "Type",
      "Injection",
      withheld_cols = character()
    ),
    "begin with a QC sample"
  )
})

test_that("clean_data errors when last sample after sort is not QC", {
  df <- data.frame(
    SampleID  = c("s1", "s2", "s3"),
    BatchID   = 1,
    Type      = c("QC", "QC", "sample"),
    Injection = c(1, 2, 3),
    met1      = c(1, 2, 3)
  )
  expect_error(
    clean_data(
      df,
      "SampleID",
      "BatchID",
      "Type",
      "Injection",
      withheld_cols = character()
    ),
    "end with a QC sample"
  )
})

test_that("withheld cols are not present in df but are tracked", {
  df <- data.frame(
    SampleID  = c("s1", "s2", "s3", "s4"),
    BatchID   = 1,
    Type      = c("QC", "QC", "sample", "QC"),
    Injection = c(1, 2, 3, 4),
    metA      = c(0, 1, 2, 3),
    keep_me   = c("a", "b", "c", "d")
  )
  out <- clean_data(df,
                    "SampleID",
                    "BatchID",
                    "Type",
                    "Injection",
                    withheld_cols = "keep_me")
  expect_false("keep_me" %in% names(out$df))
  expect_equal(out$withheld_cols, "keep_me")
})

test_that("clean_data creates a default batch column when no batch column is supplied", {
  df <- data.frame(
    SampleID = paste0("s", 1:4),
    Type = c("QC", "sample", "sample", "QC"),
    Injection = 1:4,
    met1 = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  out <- clean_data(
    df,
    sample = "SampleID",
    batch = "MissingBatch",
    class = "Type",
    order = "Injection",
    withheld_cols = character()
  )

  expect_true("batch" %in% names(out$df))
  expect_true(all(out$df$batch == "batch1"))
  expect_true("batch" %in% names(out$meta_df))
  expect_true(all(out$meta_df$batch == "batch1"))
})

test_that("clean_data removes HP rows before QC boundary checks", {
  df <- data.frame(
    SampleID = paste0("s", 1:6),
    BatchID = 1,
    Type = c("HP", "QC", "sample", "sample", "QC", "HP"),
    Injection = 1:6,
    met1 = c(99, 1, 2, 3, 4, 99),
    stringsAsFactors = FALSE
  )

  out <- clean_data(
    df,
    sample = "SampleID",
    batch = "BatchID",
    class = "Type",
    order = "Injection",
    withheld_cols = character()
  )

  expect_false(any(out$df$class == "HP"))
  expect_identical(out$df$sample, c("s2", "s3", "s4", "s5"))
  expect_identical(out$df$class[1], "QC")
  expect_identical(out$df$class[nrow(out$df)], "QC")
})

test_that("clean_data reports and repairs duplicate column names", {
  df <- data.frame(
    SampleID = paste0("s", 1:4),
    BatchID = 1,
    Type = c("QC", "sample", "sample", "QC"),
    Injection = 1:4,
    met = c(1, 2, 3, 4),
    met = c(1, 2, 3, 4),
    check.names = FALSE
  )

  out <- clean_data(
    df,
    sample = "SampleID",
    batch = "BatchID",
    class = "Type",
    order = "Injection",
    withheld_cols = character()
  )

  expect_identical(out$duplicate_col_names, "met")
  expect_true(all(c("met", "met_1") %in% names(out$df)))
})

test_that("clean_data removes entirely non-numeric metabolite columns", {
  df <- data.frame(
    SampleID = paste0("s", 1:4),
    BatchID = 1,
    Type = c("QC", "sample", "sample", "QC"),
    Injection = 1:4,
    met_numeric = c(1, 2, 3, 4),
    met_text = c("low", "medium", "high", "low"),
    stringsAsFactors = FALSE
  )

  out <- clean_data(
    df,
    sample = "SampleID",
    batch = "BatchID",
    class = "Type",
    order = "Injection",
    withheld_cols = character()
  )

  expect_false("met_text" %in% names(out$df))
  expect_identical(out$non_numeric_cols, "met_text")
  expect_identical(out$replacement_counts$metabolite, "met_numeric")
})

test_that("clean_data computes duplicate_mets correctly (nearly equal columns ignoring NAs)", {
  # Construct data so:
  # - metA and metB are equal on all rows where both non-NA (one NA mismatch allowed)
  # - metA and metC are NOT equal (different on overlapping non-NA rows)
  # Ensure QC at start/end after ordering
  df <- data.frame(
    SampleID  = paste0("s", 1:6),
    BatchID   = 1,
    Type      = c("QC", "QC", "sample", "sample", "sample", "QC"),
    Injection = 1:6,
    metA      = c(1, 2, 3, 4, NA, 6),
    metB      = c(1, 2, 3, 4, 5, 6),      # equal to metA wherever metA not NA
    metC      = c(1, 2, 30, 4, 5, 6),     # differs at Injection 3
    stringsAsFactors = FALSE
  )
  
  out <- clean_data(
    df,
    sample = "SampleID",
    batch  = "BatchID",
    class  = "Type",
    order  = "Injection",
    withheld_cols = character()
  )
  
  dm <- out$duplicate_mets
  expect_true(is.data.frame(dm))
  expect_named(dm, c("col1", "col2"))
  
  # Only metA~metB should be flagged
  expect_equal(nrow(dm), 1L)
  expect_identical(dm$col1[1], "metA")
  expect_identical(dm$col2[1], "metB")
})

test_that("clean_data duplicate_mets is empty when no equal pairs exist", {
  df <- data.frame(
    SampleID  = paste0("s", 1:5),
    BatchID   = 1,
    Type      = c("QC", "sample", "sample", "sample", "QC"),
    Injection = 1:5,
    met1      = c(1, 2, 3, 4, 5),
    met2      = c(1, 2, 3, 4, 6),  # differs at last
    met3      = c(5, 4, 3, 2, 1),
    stringsAsFactors = FALSE
  )
  
  out <- clean_data(
    df,
    sample = "SampleID",
    batch  = "BatchID",
    class  = "Type",
    order  = "Injection",
    withheld_cols = character()
  )
  
  expect_true(is.data.frame(out$duplicate_mets))
  expect_equal(nrow(out$duplicate_mets), 0L)
  expect_setequal(names(out$duplicate_mets), c("col1", "col2"))
})

test_that("clean_data separates blanks for downstream blank-threshold detection", {
  df <- data.frame(
    SampleID  = paste0("s", 1:8),
    BatchID   = 1,
    Type      = c("QC", "blank", "sample", "sample", "blank", "sample", "sample", "QC"),
    Injection = 1:8,
    # QC rows are s1 and s8
    # blanks are s2 and s5 with mean = 1 for both metabolites
    met_high  = c(10, 1,  2,  2, 1,  2,  2, 10),  # QC mean = 10, blank mean = 1 => 10 >= 3 OK
    met_low   = c( 2, 1, 10, 10, 1, 10, 10,  2),  # QC mean = 2,  blank mean = 1 => 2 < 3 FLAG
    stringsAsFactors = FALSE
  )
  
  out <- clean_data(
    df,
    sample = "SampleID",
    batch  = "BatchID",
    class  = "Type",
    order  = "Injection",
    withheld_cols = character()
  )
  
  expect_true(is.data.frame(out$blank_df))
  expect_equal(nrow(out$blank_df), 2L)
  expect_true(all(tolower(trimws(out$blank_df$class)) == "blank"))
  expect_false(any(tolower(trimws(out$df$class)) == "blank"))

  threshold <- detect_blank_threshold(
    df = out$df,
    blank_df = out$blank_df,
    metab_cols = c("met_high", "met_low")
  )

  expect_type(threshold$below_blank_threshold, "character")
  expect_true("met_low" %in% threshold$below_blank_threshold)
  expect_false("met_high" %in% threshold$below_blank_threshold)
})

test_that("clean_data returns empty blank_df when no blanks exist", {
  df <- data.frame(
    SampleID  = paste0("s", 1:5),
    BatchID   = 1,
    Type      = c("QC", "sample", "sample", "sample", "QC"),
    Injection = 1:5,
    met1      = c(1, 2, 3, 4, 5),
    met2      = c(10, 10, 10, 10, 10),
    stringsAsFactors = FALSE
  )
  
  out <- clean_data(
    df,
    sample = "SampleID",
    batch  = "BatchID",
    class  = "Type",
    order  = "Injection",
    withheld_cols = character()
  )
  
  expect_true(is.data.frame(out$blank_df))
  expect_equal(nrow(out$blank_df), 0L)
})

test_that("detect_blank_threshold does not flag metabolites when blank mean is 0 or non-finite", {
  df <- data.frame(
    SampleID  = paste0("s", 1:7),
    BatchID   = 1,
    Type      = c("QC", "blank", "sample", "sample", "blank", "sample", "QC"),
    Injection = 1:7,
    # blanks at s2 and s5 => mean = 0
    met_zero_blank = c(1, 0, 9, 9, 0, 9, 1),  # blank mean 0 => never flagged
    # blanks mean = 2, QCs are 2 and 2 => QC mean 2 < 6 => flagged
    met_flag       = c(2, 2, 100, 100, 2, 100, 2),
    stringsAsFactors = FALSE
  )
  
  out <- clean_data(
    df,
    sample = "SampleID",
    batch  = "BatchID",
    class  = "Type",
    order  = "Injection",
    withheld_cols = character()
  )

  threshold <- detect_blank_threshold(
    df = out$df,
    blank_df = out$blank_df,
    metab_cols = c("met_zero_blank", "met_flag")
  )

  expect_false("met_zero_blank" %in% threshold$below_blank_threshold)
  expect_true("met_flag" %in% threshold$below_blank_threshold)
})

test_that("repair_duplicate_column_names avoids generated-name collisions", {
  df <- data.frame(
    a = 1,
    a_1 = 2,
    a = 3,
    check.names = FALSE
  )

  out <- repair_duplicate_column_names(df)

  expect_identical(names(out), c("a", "a_1", "a_2"))
})


