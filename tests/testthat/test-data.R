test_that("embedded Hilaria jamesii data have the expected shape", {
  expect_s4_class(hija_raster, "RasterLayer")
  expect_true(is.data.frame(hija_struc))
  expect_true(ncol(hija_struc) >= 4)
  expect_true(all(c("V2", "V3") %in% names(hija_struc)))
  expect_equal(nrow(hija_struc), 16)
  expect_true(is.data.frame(hija_herb))
  expect_equal(ncol(hija_herb), 2)
})

test_that("least-cost mode reports a clear optional dependency error", {
  if (requireNamespace("gdistance", quietly = TRUE)) {
    skip("gdistance is installed")
  }

  ex_raster <- raster::aggregate(hija_raster, fact = 64)

  expect_error(
    popmaps(
      input_raster = ex_raster,
      input_locs = hija_struc,
      surface = "C",
      ncore = 1
    ),
    "gdistance"
  )
})
