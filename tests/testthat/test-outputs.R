test_that("popmaps_rast returns named terra layers aligned to the input raster", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  result <- popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    surface = "G",
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.05,
    threshold = 0,
    ncore = 1
  )

  raster_result <- popmaps_rast(result, ex_raster)

  expect_s4_class(raster_result, "SpatRaster")
  expect_equal(terra::nrow(raster_result), raster::nrow(ex_raster))
  expect_equal(terra::ncol(raster_result), raster::ncol(ex_raster))
  expect_equal(terra::nlyr(raster_result), length(result))
  expect_equal(
    names(raster_result),
    c("hard_boundary", "ancestry_probability", "axis_1", "axis_2", "axis_3")
  )
  expect_equal(
    terra::values(raster_result[[1]], mat = FALSE),
    as.vector(t(result[[1]]))
  )
})

test_that("popmaps_rast validates dimensions and custom layer names", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  result <- replicate(5, matrix(1, raster::nrow(ex_raster), raster::ncol(ex_raster)), simplify = FALSE)

  named <- popmaps_rast(
    result,
    ex_raster,
    layer_names = c("boundary", "probability", "axis one", "axis two", "axis three")
  )

  expect_equal(names(named), c("boundary", "probability", "axis.one", "axis.two", "axis.three"))

  result[[1]] <- matrix(1, nrow = 1, ncol = 1)
  expect_error(
    popmaps_rast(result, ex_raster),
    "dimensions"
  )
})

test_that("write_popmaps writes one GeoTIFF per output layer", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  result <- popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    surface = "G",
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.05,
    threshold = 0,
    ncore = 1
  )
  output_dir <- file.path(tempdir(), paste0("popmaps-output-", Sys.getpid()))

  manifest <- write_popmaps(result, ex_raster, dir = output_dir, prefix = "hija")

  expect_equal(nrow(manifest), length(result))
  expect_true(all(file.exists(manifest$path)))
  expect_true(all(grepl("\\.tif$", manifest$path)))

  written <- terra::rast(manifest$path[1])
  expect_equal(terra::nrow(written), raster::nrow(ex_raster))
  expect_equal(terra::ncol(written), raster::ncol(ex_raster))

  expect_error(
    write_popmaps(result, ex_raster, dir = output_dir, prefix = "hija"),
    "already exist"
  )
})
