
test_that("check allowed and disallowed variant strings", {

  # read test cases from the testdata directory. readr::read_csv is used
  # (rather than read.csv) as it transparently strips any byte-order mark, which
  # would otherwise corrupt the first column name
  file_path <- testthat::test_path("testdata", "variant_string_tests.csv")
  df_string <- readr::read_csv(file_path, show_col_types = FALSE)

  # guard against a silent fixture-loading failure: if the file cannot be read
  # as expected then fail loudly here rather than testing nothing
  testthat::expect_true(all(c("string", "allowed") %in% names(df_string)))
  testthat::expect_gt(nrow(df_string), 0)

  # store if each string value passes checks
  pass_checks <- rep(NA, nrow(df_string))
  for (i in 1:nrow(df_string)) {
    pass_checks[i] <- tryCatch(
      {
        # Attempt to run the function
        check_variant_string(df_string$string[i]) |>
          suppressMessages() |>
          suppressWarnings()
        TRUE  # If no error, return TRUE
      },
      error = function(e) {
        # Handle the error
        FALSE  # If there's an error, return FALSE
      }
    )
  }

  # check that this matches with what we expect
  testthat::expect_equal(pass_checks, df_string$allowed)

})
