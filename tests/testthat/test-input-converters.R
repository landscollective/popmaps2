test_that("surface_from_points builds a geographic prediction surface", {
  locs <- data.frame(
    site = c("a", "b", "c"),
    lon = c(0, 2, 4),
    lat = c(0, 1, 0),
    axis1 = c(1, 0.5, 0),
    axis2 = c(0, 0.5, 1)
  )

  surface <- surface_from_points(locs, resolution = 1, buffer = 1, surface = "G")

  expect_s3_class(surface, "popmaps_surface")
  expect_equal(surface$surface, "G")
  expect_true(is.na(surface$surface_values))
  expect_gte(terra::ncell(surface$rast), 1)
})

test_that("surface_from_points can build a constant conductance surface", {
  coords <- data.frame(x = c(0, 2), y = c(0, 2))

  surface <- surface_from_points(
    coords,
    resolution = 1,
    buffer = 0,
    surface = "C",
    surface_values = "conductance",
    values = 2
  )

  expect_equal(surface$surface, "C")
  expect_equal(surface$surface_values, "conductance")
  expect_true(all(terra::values(surface$conductance, mat = FALSE) == 2))
})

test_that("locs_from_sf converts point features to POPMAPS location tables", {
  testthat::skip_if_not_installed("sf")

  points <- sf::st_as_sf(
    data.frame(
      site = c("a", "b"),
      axis1 = c(1, 0),
      axis2 = c(0, 1),
      lon = c(-110, -111),
      lat = c(39, 40)
    ),
    coords = c("lon", "lat"),
    crs = 4326
  )

  locs <- locs_from_sf(points)

  expect_equal(names(locs), c("site", "lon", "lat", "axis1", "axis2"))
  expect_equal(locs$site, c("a", "b"))
  expect_equal(locs$axis1, c(1, 0))
  expect_equal(locs$axis2, c(0, 1))
  expect_true(all(is.finite(locs$lon)))
  expect_true(all(is.finite(locs$lat)))
})

test_that("surfaces_from_raster_stack prepares one candidate per layer", {
  rast <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(rast[[1]]) <- c(1, 2, 3, 4)
  terra::values(rast[[2]]) <- c(4, 3, 2, 1)
  names(rast) <- c("suitability", "resistance")

  surfaces <- surfaces_from_raster_stack(
    rast,
    surface = "C",
    surface_values = c("conductance", "resistance"),
    include_geographic = TRUE
  )

  expect_named(surfaces, c("geographic", "suitability", "resistance"))
  expect_equal(surfaces$geographic$surface, "G")
  expect_equal(surfaces$suitability$surface_values, "conductance")
  expect_equal(surfaces$resistance$surface_values, "resistance")
  expect_equal(surfaces$resistance$transform, "inverse_resistance")
})

test_that("surface_from_eems rasterizes coordinate-value tables", {
  values <- data.frame(
    x = c(0.5, 1.5, 0.5, 1.5),
    y = c(0.5, 0.5, 1.5, 1.5),
    migration = c(1, 2, 3, 4)
  )

  surface <- surface_from_eems(
    values,
    resolution = 1,
    buffer = 0.5,
    value_col = "migration",
    surface_values = "conductance"
  )

  expect_s3_class(surface, "popmaps_surface")
  expect_equal(surface$surface, "C")
  expect_equal(surface$surface_values, "conductance")
  expect_equal(sort(terra::values(surface$conductance, mat = FALSE)), 1:4)
})

test_that("surface_from_feems accepts FEEMS-style weight columns", {
  values <- data.frame(
    lon = c(0.5, 1.5, 0.5, 1.5),
    lat = c(0.5, 0.5, 1.5, 1.5),
    w = c(0.5, 1, 2, 4)
  )

  surface <- surface_from_feems(
    values,
    resolution = 1,
    buffer = 0.5,
    surface_values = "conductance"
  )

  expect_s3_class(surface, "popmaps_surface")
  expect_equal(surface$surface, "C")
  expect_equal(surface$surface_values, "conductance")
  expect_equal(sort(terra::values(surface$conductance, mat = FALSE)), c(0.5, 1, 2, 4))
})
