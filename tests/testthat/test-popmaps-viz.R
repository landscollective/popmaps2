popmaps_test_output <- function(input_raster, n_axes) {
  nrows <- raster::nrow(input_raster)
  ncols <- raster::ncol(input_raster)
  cell_ids <- matrix(seq_len(nrows * ncols), nrow = nrows, ncol = ncols)
  boundary <- matrix(((cell_ids - 1L) %% n_axes) + 1L, nrow = nrows, ncol = ncols)
  probability <- matrix(seq(0.1, 0.95, length.out = nrows * ncols), nrow = nrows, ncol = ncols)
  axes <- lapply(seq_len(n_axes), function(idx) {
    matrix(idx / n_axes, nrow = nrows, ncol = ncols)
  })

  c(list(boundary, probability), axes)
}

test_that("plot_popmaps draws supported modern map types", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  result <- popmaps_test_output(ex_raster, n_axes = 3)
  output_path <- file.path(tempdir(), paste0("popmaps-viz-", Sys.getpid(), ".png"))

  grDevices::png(output_path, width = 900, height = 700, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_s4_class(
    plot_popmaps(result, ex_raster, hija_struc, type = "ancestry", legend = FALSE),
    "SpatRaster"
  )
  expect_silent(plot_popmaps(result, ex_raster, hija_struc, type = "boundary", sites = "points", legend = FALSE))
  expect_silent(plot_popmaps(result, ex_raster, hija_struc, type = "axis", axis = 1, sites = "none", legend = FALSE))
  expect_silent(plot_popmaps(result, ex_raster, hija_struc, type = "axis", axis = "axis_2", sites = "none", legend = FALSE))
})

test_that("write_popmaps_plot writes a PNG and protects existing files", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  result <- popmaps_test_output(ex_raster, n_axes = 3)
  output_path <- file.path(tempdir(), paste0("popmaps-modern-map-", Sys.getpid(), ".png"))

  written <- write_popmaps_plot(
    result,
    ex_raster,
    output_path,
    input_locs = hija_struc,
    type = "boundary",
    sites = "points",
    legend = FALSE,
    width = 900,
    height = 700,
    res = 140
  )

  expect_true(file.exists(written))
  expect_gt(file.info(written)$size, 0)
  expect_error(
    write_popmaps_plot(result, ex_raster, output_path, overwrite = FALSE),
    "already exists"
  )
})

test_that("plot_popmaps validates axis and site inputs", {
  ex_raster <- raster::aggregate(hija_raster, fact = 240)
  result <- popmaps_test_output(ex_raster, n_axes = 3)
  two_axis_locs <- hija_struc[, 1:5]

  expect_error(
    plot_popmaps(result, ex_raster, type = "axis", axis = 99),
    "available ancestry-axis"
  )
  expect_error(
    plot_popmaps(result, ex_raster, two_axis_locs, sites = "points"),
    "must match"
  )
})
