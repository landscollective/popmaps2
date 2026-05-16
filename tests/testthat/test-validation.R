test_that("baseline validation matches the POPMAPS 1.03 reference", {
  validation <- validate_popmaps_baseline(quiet = TRUE)

  expect_true(isTRUE(attr(validation, "passed")))
  expect_true(all(validation$passed))
  expect_equal(validation$max_abs_diff, rep(0, nrow(validation)))
})
