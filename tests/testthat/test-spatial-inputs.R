test_that("raster inputs are normalized through terra", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  prepared <- popmaps2:::popmaps_prepare_raster(terra::rast(ex_raster))

  expect_s4_class(prepared$rast, "SpatRaster")
  expect_s4_class(prepared$raster, "RasterLayer")
  expect_equal(dim(prepared$raster), dim(ex_raster))
})

test_that("modern cell coordinates use actual raster cell centers", {
  r <- raster::raster(nrows = 3, ncols = 4, xmn = 0, xmx = 4, ymn = 0, ymx = 3)

  modern <- popmaps2:::popmaps_cell_coords(r)
  legacy <- popmaps2:::popmaps_cell_coords(r, legacy_compat = TRUE)

  expect_equal(modern, raster::xyFromCell(r, seq_len(raster::ncell(r))), ignore_attr = TRUE)
  expect_equal(modern[1, ], c(0.5, 2.5), ignore_attr = TRUE)
  expect_equal(modern[raster::ncell(r), ], c(3.5, 0.5), ignore_attr = TRUE)
  expect_false(identical(modern, legacy))
  expect_true(any(duplicated(legacy)))
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

  blank_site <- hija_struc
  blank_site$V1[1] <- ""
  expect_error(
    popmaps2:::popmaps_prepare_locations(blank_site),
    "site names"
  )

  duplicated_site <- hija_struc
  duplicated_site$V1[2] <- duplicated_site$V1[1]
  expect_error(
    popmaps2:::popmaps_prepare_locations(duplicated_site),
    "unique"
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

test_that("ancestry coefficients are checked before modeling starts", {
  zero_ancestry <- hija_struc
  zero_ancestry[1, 4:ncol(zero_ancestry)] <- 0
  expect_error(
    popmaps2:::popmaps_prepare_locations(zero_ancestry),
    "sum to a positive value"
  )

  non_normalized <- hija_struc
  non_normalized[1, 4:ncol(non_normalized)] <- non_normalized[1, 4:ncol(non_normalized)] * 0.5
  expect_warning(
    popmaps2:::popmaps_prepare_locations(non_normalized),
    "usually should sum to 1"
  )
})

test_that("empirical coordinates must fall on valid raster cells", {
  outside_locs <- hija_struc
  outside_locs$V2[1] <- raster::xmin(hija_raster) - 1

  expect_error(
    popmaps2:::popmaps_prepare_inputs(
      input_raster = hija_raster,
      input_locs = outside_locs,
      num_sites = 5,
      num_tested = 2
    ),
    "outside the raster extent"
  )

  na_raster <- terra::rast(hija_raster)
  first_cell <- terra::cellFromXY(na_raster, as.matrix(hija_struc[1, 2:3, drop = FALSE]))
  na_values <- terra::values(na_raster, mat = FALSE)
  na_values[first_cell] <- NA_real_
  terra::values(na_raster) <- na_values

  expect_error(
    popmaps2:::popmaps_prepare_inputs(
      input_raster = na_raster,
      input_locs = hija_struc,
      num_sites = 5,
      num_tested = 2
    ),
    "NA` raster cells"
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

test_that("popmaps can report progress for longer-running steps", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  messages <- character()
  result <- withCallingHandlers(
    popmaps(
      input_raster = ex_raster,
      input_locs = hija_struc,
      surface = "G",
      empirical_pt_dist = 0,
      num_sites = 5,
      num_tested = 2,
      popmod = -0.05,
      threshold = 0,
      ncore = 1,
      quiet = FALSE
    ),
    message = function(msg) {
      messages <<- c(messages, conditionMessage(msg))
      invokeRestart("muffleMessage")
    }
  )

  expect_length(result, 2 + ncol(hija_struc) - 3)
  expect_true(any(grepl("Validating raster", messages)))
  expect_true(any(grepl("Preparing geographic distances", messages)))
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

test_that("legacy compatibility remains explicit and reproducible", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)

  modern <- popmaps(
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
  legacy <- popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    surface = "G",
    empirical_pt_dist = 0,
    num_sites = 5,
    num_tested = 2,
    popmod = -0.05,
    ncore = 1,
    threshold = 0,
    legacy_compat = TRUE
  )

  expect_false(identical(modern, legacy))
  expect_equal(dim(modern[[1]]), dim(legacy[[1]]))
})
