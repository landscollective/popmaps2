test_that("embedded Hilaria jamesii data have the expected shape", {
  expect_s4_class(hija_raster, "RasterLayer")
  expect_true(is.data.frame(hija_struc))
  expect_true(ncol(hija_struc) >= 4)
  expect_true(all(c("V2", "V3") %in% names(hija_struc)))
  expect_equal(nrow(hija_struc), 16)
  expect_true(is.data.frame(hija_herb))
  expect_equal(ncol(hija_herb), 2)
})

test_that("least-cost mode runs without the optional legacy gdistance package", {
  if (requireNamespace("gdistance", quietly = TRUE)) {
    skip("gdistance is installed")
  }

  ex_raster <- raster::aggregate(hija_raster, fact = 64)

  result <- popmaps(
    input_raster = ex_raster,
    input_locs = hija_struc,
    surface = "C",
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
