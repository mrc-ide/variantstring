
test_that("check allowed and disallowed variant strings", {

  # read data.frame from inst/extdata directory
  file_path <- system.file("extdata", "variant_string_tests.csv", package = "variantstring")
  df_string <- read.csv(file_path)

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
  testthat::expect_true(all.equal.formula(df_string$allowed, pass_checks))

})
