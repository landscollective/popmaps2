test_that("raster inputs are normalized through terra", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  prepared <- popmaps2:::popmaps_prepare_raster(terra::rast(ex_raster))

  expect_s4_class(prepared$rast, "SpatRaster")
  expect_s4_class(prepared$raster, "RasterLayer")
  expect_equal(dim(prepared$raster), dim(ex_raster))
})

test_that("location inputs accept descriptive names", {
  locs <- hija_struc
  names(locs) <- c("site", "lon", "lat", "axis1", "axis2", "axis3")

  prepared <- popmaps2:::popmaps_prepare_locations(locs)

  expect_equal(names(prepared), paste0("V", seq_len(ncol(locs))))
  expect_type(prepared$V2, "double")
  expect_type(prepared$V3, "double")
  expect_true(all(prepared[, 4:6] >= 0))
})

test_that("invalid locations fail before modeling starts", {
  locs <- hija_struc
  locs$V2[1] <- "not-a-longitude"

  expect_error(
    popmaps2:::popmaps_prepare_locations(locs),
    "longitude"
  )

  expect_error(
    popmaps2:::popmaps_prepare_locations(hija_struc[, 1:3]),
    "at least four columns"
  )
})

test_that("model argument validation catches impossible site counts", {
  expect_error(
    popmaps2:::popmaps_prepare_inputs(
      input_raster = hija_raster,
      input_locs = hija_struc,
      num_sites = 99,
      num_tested = 2
    ),
    "num_sites"
  )

  expect_error(
    popmaps2:::popmaps_prepare_inputs(
      input_raster = hija_raster,
      input_locs = hija_struc,
      num_sites = 2,
      num_tested = 3
    ),
    "num_tested"
  )
})

test_that("anc_extract uses terra cell lookup and validates coordinates", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  num_axes <- ncol(hija_struc) - 3
  template <- matrix(1, nrow = raster::nrow(ex_raster), ncol = raster::ncol(ex_raster))
  pop_raster_list <- c(list(template, template), rep(list(template), num_axes))

  extracted <- anc_extract(
    pop_raster_list = pop_raster_list,
    input_raster = terra::rast(ex_raster),
    input_locs = hija_struc,
    dec_lat = 39.46522,
    dec_long = -110.9525
  )

  expect_length(extracted, num_axes)
  expect_equal(sum(extracted), 1)

  expect_error(
    anc_extract(
      pop_raster_list = pop_raster_list,
      input_raster = ex_raster,
      input_locs = hija_struc,
      dec_lat = 90,
      dec_long = 0
    ),
    "outside"
  )
})

test_that("popmaps accepts terra rasters and descriptive location columns", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  locs <- hija_struc
  names(locs) <- c("site", "lon", "lat", "axis1", "axis2", "axis3")

  result <- popmaps(
    input_raster = terra::rast(ex_raster),
    input_locs = locs,
    surface = "G",
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.05,
    threshold = 0,
    ncore = 1
  )

  expect_length(result, 2 + ncol(hija_struc) - 3)
  expect_equal(dim(result[[1]]), dim(ex_raster)[1:2])
})

test_that("popmaps preserves the legacy positional argument order", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  named <- popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    surface = "G",
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.05,
    ncore = 1,
    threshold = 0
  )
  positional <- popmaps(ex_raster, hija_struc, "G", 0, 5, 2, -0.05, 1, 0)

  expect_equal(positional, named)
})
