test_that("modern least-cost distances match gdistance on conductance rasters", {
  skip_if_not_installed("gdistance")

  conductance_rast <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  terra::values(conductance_rast) <- c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  )
  surface <- prepare_popmaps_surface(
    conductance_rast,
    surface = "C",
    surface_values = "conductance"
  )
  coords <- terra::xyFromCell(conductance_rast, c(1, 3, 5, 7, 9))

  modern <- popmaps2:::popmaps_cost_distance_matrix(surface, coords, directions = 8)

  legacy_transition <- gdistance::transition(
    raster::raster(conductance_rast),
    transitionFunction = mean,
    directions = 8
  )
  legacy_transition <- gdistance::geoCorrection(legacy_transition, type = "c", scl = TRUE)
  legacy <- as.matrix(gdistance::costDistance(legacy_transition, coords))

  expect_equal(modern, legacy, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("modern least-cost distances match gdistance with NA and zero cells", {
  skip_if_not_installed("gdistance")

  conductance_rast <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
  terra::values(conductance_rast) <- c(
    1, 0, 2, 3,
    2, NA, 2, 4,
    3, 2, 0, 5,
    4, 4, 4, 6
  )
  surface <- prepare_popmaps_surface(
    conductance_rast,
    surface = "C",
    surface_values = "conductance"
  )
  from_coords <- terra::xyFromCell(conductance_rast, c(1, 4, 16))
  to_coords <- terra::xyFromCell(conductance_rast, c(3, 9, 14))

  modern <- popmaps2:::popmaps_cost_distance_matrix(
    surface,
    from_coords = from_coords,
    to_coords = to_coords,
    directions = 8
  )

  legacy_transition <- gdistance::transition(
    raster::raster(conductance_rast),
    transitionFunction = mean,
    directions = 8
  )
  legacy_transition <- gdistance::geoCorrection(legacy_transition, type = "c", scl = TRUE)
  legacy <- gdistance::costDistance(legacy_transition, from_coords, to_coords)

  expect_equal(modern, legacy, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("least-cost cell distance helper returns one row per raster cell", {
  conductance_rast <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  terra::values(conductance_rast) <- c(
    1, 1, 1,
    1, NA, 1,
    1, 1, 1
  )
  surface <- prepare_popmaps_surface(
    conductance_rast,
    surface = "C",
    surface_values = "conductance"
  )
  from_coords <- terra::xyFromCell(conductance_rast, c(1, 9))

  distances <- popmaps2:::popmaps_cost_distance_to_cells(surface, from_coords, directions = 4)

  expect_equal(dim(distances), c(9, 2))
  expect_true(all(is.na(distances[5, ])))
  expect_equal(distances[1, 1], 0)
  expect_equal(distances[9, 2], 0)
  expect_true(is.finite(distances[3, 1]))
})

test_that("least-cost helper honors resistance-to-conductance surfaces", {
  resistance_rast <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2)
  terra::values(resistance_rast) <- c(1, 2, 4, 1, 2, 4)
  surface <- prepare_popmaps_surface(
    resistance_rast,
    surface = "C",
    surface_values = "resistance",
    resistance_epsilon = 0.5
  )
  coords <- terra::xyFromCell(resistance_rast, c(1, 3))

  distances <- popmaps2:::popmaps_cost_distance_matrix(surface, coords, directions = 4)

  expect_equal(distances[1, 1], 0)
  expect_equal(distances[2, 2], 0)
  expect_true(distances[1, 2] > 0)
})

test_that("least-cost helper rejects geographic surfaces and non-traversable coordinates", {
  conductance_rast <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(conductance_rast) <- c(1, NA, 1, 1)
  geographic <- prepare_popmaps_surface(conductance_rast, surface = "G")
  conductance <- prepare_popmaps_surface(conductance_rast, surface = "C")

  expect_error(
    popmaps2:::popmaps_cost_distance_graph(geographic),
    "surface = \"C\""
  )
  expect_error(
    popmaps2:::popmaps_cost_distance_matrix(
      conductance,
      terra::xyFromCell(conductance_rast, 2)
    ),
    "traversable"
  )
})
