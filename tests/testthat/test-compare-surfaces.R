surface_comparison_fixture <- function() {
  surface_rast <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
  terra::values(surface_rast) <- c(
    1, 1, 1, 1,
    1, 2, 2, 1,
    1, 2, 2, 1,
    1, 1, 1, 1
  )
  locs <- data.frame(
    site = paste0("s", 1:4),
    lon = c(0.5, 3.5, 0.5, 3.5),
    lat = c(0.5, 0.5, 3.5, 3.5),
    axis1 = c(0.9, 0.8, 0.2, 0.1),
    axis2 = c(0.1, 0.2, 0.8, 0.9)
  )

  list(raster = surface_rast, locs = locs)
}

test_that("compare_popmaps_surfaces ranks prepared and specified surfaces", {
  fixture <- surface_comparison_fixture()
  surfaces <- list(
    geographic = prepare_popmaps_surface(fixture$raster, surface = "G"),
    suitability = list(
      input_raster = fixture$raster,
      surface = "C",
      surface_values = "suitability"
    )
  )

  comparison <- compare_popmaps_surfaces(
    input_locs = fixture$locs,
    surfaces = surfaces,
    empirical_pt_dist = c(0, 1),
    num_sites = 3,
    num_tested = 2,
    popmod = c(-0.1, -0.2),
    quiet = TRUE
  )

  expect_s3_class(comparison, "popmaps_surface_comparison")
  expect_named(
    comparison,
    c(
      "summary", "best", "near_best", "support", "tunings",
      "grids", "primary_metric", "validation", "surface_grid", "near_best_tolerance",
      "spatial_block_seed", "call"
    )
  )
  expect_equal(nrow(comparison$summary), 2)
  expect_equal(names(comparison$tunings), c("geographic", "suitability"))
  expect_equal(sort(comparison$summary$surface_name), c("geographic", "suitability"))
  expect_equal(comparison$summary$rank, 1:2)
  expect_true(all(is.finite(comparison$summary$score)))
  expect_true(all(comparison$summary$failed_folds == 0))
  expect_equal(nrow(comparison$best), 1)
  expect_true(comparison$best$surface_name %in% comparison$summary$surface_name)
  expect_true(comparison$support$n_near_best >= 1)
})

test_that("compare_popmaps_surfaces can suggest surface-specific distance grids", {
  fixture <- surface_comparison_fixture()
  surfaces <- list(
    geographic = prepare_popmaps_surface(fixture$raster, surface = "G"),
    suitability = prepare_popmaps_surface(
      fixture$raster,
      surface = "C",
      surface_values = "suitability"
    )
  )

  comparison <- compare_popmaps_surfaces(
    input_locs = fixture$locs,
    surfaces = surfaces,
    num_sites = 3,
    num_tested = 2,
    validation = "loo",
    quiet = TRUE
  )

  expect_equal(comparison$surface_grid, "surface_specific")
  expect_equal(comparison$grids$geographic$distance_units, "km")
  expect_equal(comparison$grids$suitability$distance_units, "cost_distance")
  expect_equal(
    sort(unique(comparison$tunings$geographic$results$popmod)),
    comparison$grids$geographic$popmod
  )
  expect_equal(
    sort(unique(comparison$tunings$suitability$results$popmod)),
    comparison$grids$suitability$popmod
  )
})

test_that("compare_popmaps_surfaces uses matched repeated spatial blocks", {
  fixture <- surface_comparison_fixture()
  surfaces <- list(
    geographic = prepare_popmaps_surface(fixture$raster, surface = "G"),
    suitability = prepare_popmaps_surface(
      fixture$raster,
      surface = "C",
      surface_values = "conductance"
    )
  )

  comparison <- compare_popmaps_surfaces(
    input_locs = fixture$locs,
    surfaces = surfaces,
    validation = "spatial_block",
    n_blocks = 2,
    spatial_block_repeats = 2,
    empirical_pt_dist = 0,
    num_sites = 2,
    num_tested = 2,
    popmod = -0.1,
    quiet = TRUE
  )

  expect_equal(comparison$spatial_block_seed, 1)
  expect_equal(
    comparison$tunings$geographic$folds$block_id,
    comparison$tunings$suitability$folds$block_id
  )
  expect_equal(
    comparison$tunings$geographic$folds$site_index,
    comparison$tunings$suitability$folds$site_index
  )
  expect_true(all(comparison$summary$n_validation_repeats == 2))
})

test_that("compare_popmaps_surfaces validates surface inputs", {
  fixture <- surface_comparison_fixture()

  expect_error(
    compare_popmaps_surfaces(
      input_locs = fixture$locs,
      surfaces = list(one = prepare_popmaps_surface(fixture$raster, surface = "G")),
      empirical_pt_dist = 0,
      num_sites = 3,
      num_tested = 2,
      popmod = -0.1
    ),
    "at least two"
  )

  expect_error(
    compare_popmaps_surfaces(
      input_locs = fixture$locs,
      surfaces = list(a = list(surface = "G"), b = list(surface = "C")),
      empirical_pt_dist = 0,
      num_sites = 3,
      num_tested = 2,
      popmod = -0.1
    ),
    "input_raster"
  )

  expect_error(
    compare_popmaps_surfaces(
      input_locs = fixture$locs,
      surfaces = stats::setNames(
        list(
          prepare_popmaps_surface(fixture$raster, surface = "G"),
          prepare_popmaps_surface(fixture$raster, surface = "C")
        ),
        c("same", "same")
      ),
      empirical_pt_dist = 0,
      num_sites = 3,
      num_tested = 2,
      popmod = -0.1
    ),
    "unique"
  )
})
