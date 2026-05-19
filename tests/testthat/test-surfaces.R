test_that("prepare_popmaps_surface records geographic surfaces", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  prepared <- prepare_popmaps_surface(ex_raster, surface = "G")

  expect_s3_class(prepared, "popmaps_surface")
  expect_equal(prepared$surface, "G")
  expect_true(is.na(prepared$surface_values))
  expect_null(prepared$conductance)
  expect_s4_class(prepared$rast, "SpatRaster")
})

test_that("suitability and conductance surfaces use values directly", {
  surface_rast <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(surface_rast) <- c(0, 0.25, 0.75, NA)

  prepared <- prepare_popmaps_surface(
    surface_rast,
    surface = "C",
    surface_values = "suitability"
  )

  expect_equal(prepared$surface_values, "suitability")
  expect_equal(prepared$transform, "identity")
  expect_equal(
    terra::values(prepared$conductance, mat = FALSE),
    terra::values(surface_rast, mat = FALSE)
  )
})

test_that("resistance surfaces are converted to conductance", {
  resistance_rast <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(resistance_rast) <- c(0, 1, 3, NA)

  prepared <- prepare_popmaps_surface(
    resistance_rast,
    surface = "C",
    surface_values = "resistance",
    resistance_epsilon = 0.5
  )

  expect_equal(prepared$surface_values, "resistance")
  expect_equal(prepared$transform, "inverse_resistance")
  expect_equal(
    terra::values(prepared$conductance, mat = FALSE),
    c(1 / 0.5, 1 / 1.5, 1 / 3.5, NA)
  )
})

test_that("conductance rescaling is explicit", {
  surface_rast <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(surface_rast) <- c(2, 4, 6, NA)

  prepared <- prepare_popmaps_surface(
    surface_rast,
    surface = "C",
    surface_values = "conductance",
    rescale_conductance = TRUE
  )

  expect_equal(prepared$transform, "identity+rescaled_0_1")
  expect_equal(
    terra::values(prepared$conductance, mat = FALSE),
    c(0, 0.5, 1, NA)
  )
})

test_that("mask and barrier rasters are validated separately from conductance", {
  surface_rast <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(surface_rast) <- c(1, 2, 3, 4)

  mask_rast <- surface_rast
  terra::values(mask_rast) <- c(1, 0, NA, 1)

  barrier_rast <- surface_rast
  terra::values(barrier_rast) <- c(0, 1, 0, NA)

  prepared <- prepare_popmaps_surface(
    surface_rast,
    surface = "C",
    surface_values = "conductance",
    mask = mask_rast,
    barrier = barrier_rast
  )

  expect_equal(terra::values(prepared$mask, mat = FALSE), c(1, 0, 0, 1))
  expect_equal(terra::values(prepared$barrier, mat = FALSE), c(0, 1, 0, 0))
  expect_true(is.na(terra::values(prepared$conductance, mat = FALSE)[2]))
})

test_that("optional mask and barrier rasters must match surface geometry", {
  surface_rast <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(surface_rast) <- c(1, 2, 3, 4)

  wrong_geometry <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  terra::values(wrong_geometry) <- seq_len(9)

  expect_error(
    prepare_popmaps_surface(
      surface_rast,
      surface = "C",
      surface_values = "conductance",
      mask = wrong_geometry
    ),
    "same geometry"
  )
})

test_that("invalid conductance values fail early", {
  surface_rast <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(surface_rast) <- c(1, -1, 2, NA)

  expect_error(
    prepare_popmaps_surface(
      surface_rast,
      surface = "C",
      surface_values = "conductance"
    ),
    "non-negative"
  )

  terra::values(surface_rast) <- c(0, 0, 0, NA)
  expect_error(
    prepare_popmaps_surface(
      surface_rast,
      surface = "C",
      surface_values = "conductance",
      rescale_conductance = TRUE
    ),
    "positive conductance"
  )
})
