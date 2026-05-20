gdistance_popmaps_c_reference <- function(raster_surface,
                                          species_data,
                                          empirical_pt_dist,
                                          num_sites,
                                          num_tested,
                                          popmod,
                                          threshold,
                                          dist_prob_func,
                                          legacy_compat = FALSE) {
  nrows <- raster_surface@nrows
  ncols <- raster_surface@ncols
  ymax <- raster_surface@extent@ymax
  xmin <- raster_surface@extent@xmin
  cell_size <- raster::res(raster_surface)[1]
  num_axes <- ncol(species_data) - 3
  ancestry <- as.matrix(species_data[, seq.int(4, ncol(species_data)), drop = FALSE])
  sampling_loc_coords <- as.matrix(species_data[, 2:3, drop = FALSE])
  coords <- popmaps2:::popmaps_cell_coords(
    raster_surface,
    legacy_compat = legacy_compat
  )

  raster_values <- raster::extract(raster_surface, coords)
  geographic_cell_distances <- popmaps2:::popmaps_cell_site_distances(coords, sampling_loc_coords)
  species_surface <- gdistance::transition(raster_surface, transitionFunction = mean, directions = 8)
  species_surface <- gdistance::geoCorrection(species_surface, type = "c", scl = TRUE)
  empirical_distances <- as.matrix(gdistance::costDistance(species_surface, sampling_loc_coords))
  hard_boundaries <- max.col(ancestry, ties.method = "first")
  result <- matrix(NA_real_, nrow = nrow(coords), ncol = num_axes + 2)

  for (cell_idx in seq_len(nrow(coords))) {
    if (is.na(raster_values[cell_idx])) {
      next
    }

    if (isTRUE(legacy_compat)) {
      geographic_order <- order(geographic_cell_distances[cell_idx, ])[seq_len(num_sites)]
      temp_data <- rbind(coords[cell_idx, , drop = FALSE], sampling_loc_coords[geographic_order, , drop = FALSE])
      cost_values <- as.numeric(gdistance::costDistance(species_surface, temp_data))[seq_len(num_sites)]
      cell_order <- geographic_order[order(cost_values, na.last = TRUE)]
      cost_by_site <- rep(NA_real_, nrow(sampling_loc_coords))
      cost_by_site[geographic_order] <- cost_values
    } else {
      temp_data <- rbind(coords[cell_idx, , drop = FALSE], sampling_loc_coords)
      cost_values <- as.numeric(gdistance::costDistance(species_surface, temp_data))[seq_len(nrow(sampling_loc_coords))]
      reachable_sites <- which(is.finite(cost_values))
      if (length(reachable_sites) < num_sites) {
        next
      }
      cell_order <- reachable_sites[order(cost_values[reachable_sites])][seq_len(num_sites)]
      cost_by_site <- cost_values
    }
    result[cell_idx, 1] <- hard_boundaries[cell_order[1]]

    if (raster_values[cell_idx] < threshold) {
      next
    }

    selected_sites <- popmaps2:::popmaps_select_empirical_sites(
      candidate_sites = cell_order,
      empirical_distances = empirical_distances,
      empirical_pt_dist = empirical_pt_dist,
      num_tested = num_tested
    )

    site_weights <- vapply(
      cost_by_site[selected_sites],
      function(distance) dist_prob_func(popmod, distance),
      numeric(1)
    )
    cell_prob <- colSums(ancestry[selected_sites, , drop = FALSE] * (site_weights / num_tested))
    max_avg <- max(cell_prob) / sum(cell_prob)
    if (is.na(max_avg)) {
      max_avg <- 0
    }
    if (max_avg < 0.5) {
      max_avg <- 0
    } else {
      max_avg <- (max_avg - 0.5) * 2
    }

    result[cell_idx, 2] <- max_avg
    result[cell_idx, seq.int(3, num_axes + 2)] <- cell_prob
  }

  output <- list(
    matrix(result[, 1], nrow = nrows, ncol = ncols, byrow = TRUE),
    matrix(result[, 2], nrow = nrows, ncol = ncols, byrow = TRUE)
  )
  axes <- vector("list", num_axes)
  names(axes) <- paste0("axis", seq_len(num_axes))
  for (axis_idx in seq_len(num_axes)) {
    axes[[axis_idx]] <- matrix(result[, axis_idx + 2], nrow = nrows, ncol = ncols, byrow = TRUE)
  }
  append(output, axes)
}

test_that("popmaps C surfaces match gdistance least-cost math on small rasters", {
  skip_if_not_installed("gdistance")

  surface_rast <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
  terra::values(surface_rast) <- c(
    1, 2, 3, 4,
    2, 5, 6, 3,
    3, 6, 5, 2,
    4, 3, 2, 1
  )
  locs <- data.frame(
    site = paste0("s", 1:4),
    lon = c(0.5, 3.5, 0.5, 3.5),
    lat = c(0.5, 0.5, 3.5, 3.5),
    axis1 = c(0.9, 0.8, 0.2, 0.1),
    axis2 = c(0.1, 0.2, 0.8, 0.9)
  )
  raster_surface <- raster::raster(surface_rast)

  modern <- popmaps(
    input_raster = raster_surface,
    input_locs = locs,
    surface = "C",
    surface_values = "suitability",
    empirical_pt_dist = 0,
    num_sites = 3,
    num_tested = 2,
    popmod = -0.1,
    threshold = 0,
    ncore = 1
  )
  legacy <- gdistance_popmaps_c_reference(
    raster_surface = raster_surface,
    species_data = popmaps2:::popmaps_prepare_locations(locs),
    empirical_pt_dist = 0,
    num_sites = 3,
    num_tested = 2,
    popmod = -0.1,
    threshold = 0,
    dist_prob_func = function(popmod_temp, distance) exp(popmod_temp * distance)
  )

  for (idx in seq_along(modern)) {
    expect_equal(modern[[idx]], legacy[[idx]], tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("popmaps C surfaces preserve NA cells and threshold hard boundaries", {
  surface_rast <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  terra::values(surface_rast) <- c(
    1, NA, 0.1,
    1, 1, 1,
    1, 1, 1
  )
  locs <- data.frame(
    site = paste0("s", 1:4),
    lon = c(0.5, 2.5, 0.5, 2.5),
    lat = c(0.5, 0.5, 2.5, 2.5),
    axis1 = c(0.9, 0.8, 0.2, 0.1),
    axis2 = c(0.1, 0.2, 0.8, 0.9)
  )

  result <- popmaps(
    input_raster = surface_rast,
    input_locs = locs,
    surface = "C",
    empirical_pt_dist = 0,
    num_sites = 3,
    num_tested = 2,
    popmod = -0.1,
    threshold = 0.5,
    ncore = 1
  )

  expect_true(all(is.na(vapply(result, function(layer) layer[1, 2], numeric(1)))))
  expect_false(is.na(result[[1]][1, 3]))
  expect_true(is.na(result[[2]][1, 3]))
  expect_true(is.na(result[[3]][1, 3]))
  expect_true(is.na(result[[4]][1, 3]))
})

test_that("popmaps C surfaces accept resistance rasters", {
  resistance_rast <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  terra::values(resistance_rast) <- c(
    1, 2, 4,
    1, 2, 4,
    1, 2, 4
  )
  locs <- data.frame(
    site = paste0("s", 1:4),
    lon = c(0.5, 2.5, 0.5, 2.5),
    lat = c(0.5, 0.5, 2.5, 2.5),
    axis1 = c(0.9, 0.8, 0.2, 0.1),
    axis2 = c(0.1, 0.2, 0.8, 0.9)
  )

  result <- popmaps(
    input_raster = resistance_rast,
    input_locs = locs,
    surface = "C",
    surface_values = "resistance",
    empirical_pt_dist = 0,
    num_sites = 3,
    num_tested = 2,
    popmod = -0.1,
    threshold = 0,
    ncore = 1,
    resistance_epsilon = 0.5
  )

  expect_length(result, 4)
  expect_equal(dim(result[[1]]), c(3, 3))
  expect_true(any(is.finite(result[[2]])))
})
