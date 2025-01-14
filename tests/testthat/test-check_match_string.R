
test_that("check for correct matches and mismatches of variant strings and position strings", {

  # read data.frame from inst/extdata directory
  file_path <- system.file("extdata", "match_string_tests.csv", package = "variantstring")
  df_string <- read.csv(file_path)

  # subset to variant strings
  df_variant <- df_string |>
    dplyr::filter(position_string == FALSE) |>
    dplyr::mutate(confirm_match = NA,
                  confirm_ambiguous = NA)

  # store results of each comparison
  for (i in 1:nrow(df_variant)) {
    x <- compare_variant_string(target_string = df_variant$target_string[i],
                                comparison_strings = df_variant$comparison_string[i])
    df_variant$confirm_match[i] <- x$match
    df_variant$confirm_ambiguous[i] <- x$ambiguous
  }

  # check that this matches with what we expect
  all((df_variant$expect_match == df_variant$confirm_match) &
    (df_variant$expect_ambiguous == df_variant$confirm_ambiguous)) |>
    testthat::expect_true()
})
