test_that("tune_popmaps returns parameter summaries and fold diagnostics", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  tuning <- tune_popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    empirical_pt_dist = c(0, 5),
    num_sites = c(5, 6),
    num_tested = c(2, 3),
    popmod = c(-0.01, -0.05),
    quiet = TRUE
  )

  expect_s3_class(tuning, "popmaps_tuning")
  expect_named(tuning, c("results", "folds", "best", "primary_metric", "call"))
  expect_equal(nrow(tuning$results), 16)
  expect_equal(nrow(tuning$folds), 16 * nrow(hija_struc))
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

test_that("suggest_tuning_grid builds interpretable defaults from site distances", {
  grid <- suggest_tuning_grid(hija_struc)

  expect_s3_class(grid, "popmaps_tuning_grid")
  expect_true(all(c("num_sites", "num_tested", "popmod", "empirical_pt_dist") %in% names(grid)))
  expect_true(all(grid$num_tested <= min(grid$num_sites)))
  expect_true(all(grid$popmod < 0))
  expect_true(all(grid$empirical_pt_dist >= 0))
  expect_equal(grid$empirical_pt_dist[1], 0)
  expect_true(is.finite(grid$reference_distance))
  expect_true(grid$reference_distance > 0)
})

test_that("adaptive_tune_popmaps samples and refines parameter space reproducibly", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  tuning <- adaptive_tune_popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    num_sites = c(5, 8),
    num_tested = c(2, 3),
    popmod = c(-0.02, -0.01),
    empirical_pt_dist = c(0, 5),
    n_initial = 6,
    n_refine = 4,
    method = "latin_hypercube",
    seed = 10,
    quiet = TRUE
  )

  expect_s3_class(tuning, "popmaps_adaptive_tuning")
  expect_s3_class(tuning, "popmaps_tuning")
  expect_equal(tuning$search$n_initial, 6)
  expect_equal(tuning$search$n_refine, 4)
  expect_true(nrow(tuning$results) <= 10)
  expect_true(nrow(tuning$results) >= 6)
  expect_equal(nrow(tuning$folds), nrow(tuning$results) * nrow(hija_struc))
  expect_true(all(tuning$results$num_tested <= tuning$results$num_sites))
  expect_equal(nrow(tuning$best), 1)
  expect_true(all(is.finite(tuning$results$rmse)))
})

test_that("best-row selection prefers complete validation over partial low error", {
  results <- data.frame(
    combo_id = c(1, 2),
    num_sites = c(5, 5),
    num_tested = c(2, 2),
    popmod = c(-0.1, -0.1),
    empirical_pt_dist = c(100, 0),
    n_folds = c(16, 16),
    n_scored = c(4, 16),
    failed_folds = c(12, 0),
    mae = c(0, 0.1),
    rmse = c(0, 0.1),
    hellinger = c(0, 0.1),
    dominant_accuracy = c(1, 0.8),
    dominant_probability = c(1, 0.8)
  )

  best <- popmaps2:::popmaps_tuning_best_row(results, "rmse")

  expect_equal(best$combo_id, 2)
})
