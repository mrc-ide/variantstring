test_that("check unique speed ups", {

  # read data.frame from inst/extdata directory
  file_path <- system.file("extdata", "variant_string_tests.csv", package = "variantstring")
  df_string <- read.csv(file_path)

  # make suitable long test case
  examples <- df_string$string[df_string$allowed]
  examples <- rep(examples, 100)

  # New version - should use tictoc::tic() but would need a new Suggests
  t1 <- Sys.time()
  out <- order_variant_string(examples)
  t2 <- Sys.time()
  t_now <- t2 - t1

  # Old order_variant_string
  old_order_variant_string <- function(x) {

    # checks
    check_variant_string(x)

    mapply(function(y) {
      dplyr::arrange(y, gene, pos, aa)
    }, variant_to_long(x), SIMPLIFY = FALSE) |>
      long_to_variant()

  }

  t3 <- Sys.time()
  out2 <- old_order_variant_string(examples)
  t4 <- Sys.time()
  t_old <- t4 - t3

  # check it's quicker and the same
  testthat::expect_true(t_old > t_now)
  testthat::expect_equal(out > out2)

})
