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
  expect_named(tuning, c("results", "folds", "best", "primary_metric", "validation", "call"))
  expect_equal(nrow(tuning$results), 16)
  expect_equal(nrow(tuning$folds), 16 * nrow(hija_struc))
  expect_equal(nrow(tuning$best), 1)
  expect_equal(tuning$validation, "loo")
  expect_true(all(tuning$results$n_scored == nrow(hija_struc)))
  expect_true(all(tuning$results$failed_folds == 0))
  expect_true(all(tuning$results$n_validation_repeats == 1))
  expect_true(all(tuning$results$n_validation_folds == nrow(hija_struc)))
  expect_true(all(tuning$results$half_distance_km > 0))
  expect_true(all(tuning$results$ten_pct_distance_km > tuning$results$half_distance_km))
  expect_true(all(is.finite(tuning$results$rmse)))
  expect_true(all(is.finite(tuning$folds$predicted_axis_1)))
  expect_true(all(is.finite(tuning$folds$observed_axis_1)))
  expect_equal(
    tuning$results$half_distance_km,
    log(0.5) / tuning$results$popmod
  )
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

test_that("tune_popmaps supports spatial-block validation", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  tuning <- tune_popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    validation = "spatial_block",
    n_blocks = 4,
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.01,
    quiet = TRUE
  )

  expect_s3_class(tuning, "popmaps_tuning")
  expect_equal(tuning$validation, "spatial_block")
  expect_equal(nrow(tuning$results), 1)
  expect_equal(nrow(tuning$folds), nrow(hija_struc))
  expect_true(length(unique(tuning$folds$fold_id)) >= 2)
  expect_true(all(tuning$folds$validation == "spatial_block"))
  expect_true(all(is.finite(tuning$folds$n_training)))
  expect_true(all(tuning$results$n_validation_repeats == 1))
  expect_true(all(tuning$results$n_validation_folds >= 2))
})

test_that("tune_popmaps supports repeated spatial-block validation", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  tuning <- tune_popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    validation = "spatial_block",
    n_blocks = 4,
    spatial_block_repeats = 3,
    spatial_block_seed = 99,
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.01,
    quiet = TRUE
  )
  tuning_again <- tune_popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    validation = "spatial_block",
    n_blocks = 4,
    spatial_block_repeats = 3,
    spatial_block_seed = 99,
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.01,
    quiet = TRUE
  )

  expect_equal(nrow(tuning$results), 1)
  expect_equal(nrow(tuning$folds), 3 * nrow(hija_struc))
  expect_equal(sort(unique(tuning$folds$repeat_id)), 1:3)
  expect_equal(tuning$results$n_validation_repeats, 3)
  expect_true(is.finite(tuning$results$rmse_repeat_sd))
  expect_equal(tuning$folds$block_id, tuning_again$folds$block_id)
  expect_equal(tuning$folds$rmse, tuning_again$folds$rmse)
})

test_that("tune_popmaps accepts supplied spatial block assignments", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  blocks <- rep(c("west", "east"), length.out = nrow(hija_struc))

  tuning <- tune_popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    validation = "spatial_block",
    block_assignments = blocks,
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.01,
    quiet = TRUE
  )

  expect_equal(sort(unique(tuning$folds$block_id)), c("east", "west"))
  expect_error(
    tune_popmaps(
      input_raster = ex_raster,
      input_locs = hija_struc,
      validation = "spatial_block",
      block_assignments = blocks,
      spatial_block_repeats = 2,
      empirical_pt_dist = 0,
      num_sites = 5,
      num_tested = 2,
      popmod = -0.01,
      quiet = TRUE
    ),
    "block_assignments = NULL"
  )
  expect_error(
    tune_popmaps(
      input_raster = ex_raster,
      input_locs = hija_struc,
      validation = "spatial_block",
      block_assignments = blocks[-1],
      empirical_pt_dist = 0,
      num_sites = 5,
      num_tested = 2,
      popmod = -0.01,
      quiet = TRUE
    ),
    "one value per empirical site"
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

test_that("diagnose_tuning summarizes tuning support", {
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
  diagnostics <- diagnose_tuning(tuning)

  expect_s3_class(diagnostics, "popmaps_tuning_diagnostics")
  expect_named(diagnostics, c("overview", "near_best", "parameter_ranges", "parameter_effects"))
  expect_equal(nrow(diagnostics$overview), 1)
  expect_equal(diagnostics$overview$primary_metric, "rmse")
  expect_equal(diagnostics$overview$metric_goal, "minimize")
  expect_true(diagnostics$overview$n_near_best >= 1)
  expect_true(tuning$best$combo_id %in% diagnostics$near_best$combo_id)
  expect_true(all(c("num_sites", "num_tested", "popmod", "empirical_pt_dist") %in% diagnostics$parameter_ranges$parameter))
  expect_true(all(c("parameter", "value", "mean_score", "loss_from_best", "rank") %in% names(diagnostics$parameter_effects)))
  expect_true(all(diagnostics$parameter_effects$loss_from_best >= -sqrt(.Machine$double.eps), na.rm = TRUE))

  expect_error(
    diagnose_tuning(tuning, near_best_tolerance = -0.1),
    "non-negative"
  )
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
