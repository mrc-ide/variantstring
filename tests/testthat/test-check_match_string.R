
test_that("check for correct matches and mismatches of variant strings and position strings", {

  # read test cases from the testdata directory. readr::read_csv is used
  # (rather than read.csv) as it transparently strips any byte-order mark, which
  # would otherwise corrupt the first column name
  file_path <- testthat::test_path("testdata", "match_string_tests.csv")
  df_string <- readr::read_csv(file_path, show_col_types = FALSE)

  # guard against a silent fixture-loading failure: if the file cannot be read
  # as expected then fail loudly here rather than testing nothing
  testthat::expect_true(all(c("target_string", "comparison_string", "position_string",
                              "expect_match", "expect_ambiguous", "expect_prop") %in% names(df_string)))
  testthat::expect_gt(nrow(df_string), 0)

  # subset to variant strings
  df_variant <- df_string |>
    dplyr::filter(position_string == FALSE) |>
    dplyr::mutate(confirm_match = NA,
                  confirm_ambiguous = NA,
                  confirm_prop = NA)

  # store results of each comparison
  for (i in 1:nrow(df_variant)) {
    x <- compare_variant_string(target_string = df_variant$target_string[i],
                                comparison_strings = df_variant$comparison_string[i])
    df_variant$confirm_match[i] <- x$match
    df_variant$confirm_ambiguous[i] <- x$ambiguous
    df_variant$confirm_prop[i] <- x$prop
  }

  # check that this matches with what we expect
  match_correct <- isTRUE(all.equal(df_variant$expect_match, df_variant$confirm_match))
  ambiguous_correct <- isTRUE(all.equal(df_variant$expect_ambiguous, df_variant$confirm_ambiguous))
  prop_correct <- isTRUE(all.equal(df_variant$expect_prop, df_variant$confirm_prop))

  testthat::expect_true(match_correct & ambiguous_correct & prop_correct)
})
