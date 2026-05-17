test_that("tune_popmaps returns parameter summaries and fold diagnostics", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  tuning <- tune_popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    empirical_pt_dist = 0,
    num_sites = c(5, 6),
    num_tested = c(2, 3),
    popmod = c(-0.01, -0.05),
    quiet = TRUE
  )

  expect_s3_class(tuning, "popmaps_tuning")
  expect_named(tuning, c("results", "folds", "best", "primary_metric", "call"))
  expect_equal(nrow(tuning$results), 8)
  expect_equal(nrow(tuning$folds), 8 * nrow(hija_struc))
  expect_equal(nrow(tuning$best), 1)
  expect_true(all(tuning$results$n_scored == nrow(hija_struc)))
  expect_true(all(tuning$results$failed_folds == 0))
  expect_true(all(is.finite(tuning$results$rmse)))
  expect_true(all(is.finite(tuning$folds$predicted_axis_1)))
  expect_true(all(is.finite(tuning$folds$observed_axis_1)))
})

test_that("tune_popmaps validates unsupported and impossible tuning requests", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  expect_error(
    tune_popmaps(
      input_raster = ex_raster,
      input_locs = hija_struc,
      surface = "C",
      quiet = TRUE
    ),
    "surface = 'G'"
  )

  expect_error(
    tune_popmaps(
      input_raster = ex_raster,
      input_locs = hija_struc,
      num_sites = 2,
      num_tested = 3,
      quiet = TRUE
    ),
    "num_tested <= num_sites"
  )
})
