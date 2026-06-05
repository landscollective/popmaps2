popmaps_geographic_surface <- function(raster_surface,
                                       species_data,
                                       empirical_pt_dist,
                                       num_sites,
                                       num_tested,
                                       popmod,
                                       threshold,
                                       dist_prob_func,
                                       legacy_compat = FALSE,
                                       quiet = TRUE) {
  nrows <- raster_surface@nrows
  ncols <- raster_surface@ncols

  num_axes <- ncol(species_data) - 3
  ancestry <- as.matrix(species_data[, seq.int(4, ncol(species_data)), drop = FALSE])
  sampling_loc_coords <- as.matrix(species_data[, 2:3, drop = FALSE])

  coords <- popmaps_cell_coords(raster_surface, legacy_compat = legacy_compat)
  raster_values <- raster::extract(raster_surface, coords)
  popmaps_inform(
    "Preparing geographic distances for ",
    nrow(coords),
    " raster cells and ",
    nrow(sampling_loc_coords),
    " empirical sites.",
    quiet = quiet
  )
  cell_distances <- popmaps_cell_site_distances(coords, sampling_loc_coords)
  empirical_distances <- popmaps_empirical_site_distances(sampling_loc_coords)
  hard_boundaries <- max.col(ancestry, ties.method = "first")

  result <- matrix(NA_real_, nrow = nrow(coords), ncol = num_axes + 2)

  for (cell_idx in seq_len(nrow(coords))) {
    popmaps_progress_message(cell_idx, nrow(coords), "Estimated ancestry for raster cell", quiet = quiet)
    if (is.na(raster_values[cell_idx])) {
      next
    }

    cell_order <- order(cell_distances[cell_idx, ])[seq_len(num_sites)]
    nearest_site <- cell_order[1]
    result[cell_idx, 1] <- hard_boundaries[nearest_site]

    if (raster_values[cell_idx] < threshold) {
      next
    }

    selected_sites <- popmaps_select_empirical_sites(
      candidate_sites = cell_order,
      empirical_distances = empirical_distances,
      empirical_pt_dist = empirical_pt_dist,
      num_tested = num_tested
    )

    site_weights <- vapply(
      cell_distances[cell_idx, selected_sites],
      function(distance) dist_prob_func(popmod, distance),
      numeric(1)
    )

    weighted_ancestry <- colSums(ancestry[selected_sites, , drop = FALSE] * (site_weights / num_tested))
    dominant_share <- max(weighted_ancestry) / sum(weighted_ancestry)

    if (is.na(dominant_share)) {
      dominant_share <- 0
    }

    # Convert dominant ancestry share to the legacy POPMAPS confidence scale.
    # This reports assignment dominance, not posterior probability.
    if (dominant_share < 0.5) {
      dominant_confidence <- 0
    } else {
      dominant_confidence <- (dominant_share - 0.5) * 2
    }

    result[cell_idx, 2] <- dominant_confidence
    result[cell_idx, seq.int(3, num_axes + 2)] <- weighted_ancestry
  }

  output <- list(
    matrix(result[, 1], nrow = nrows, ncol = ncols, byrow = TRUE),
    matrix(result[, 2], nrow = nrows, ncol = ncols, byrow = TRUE)
  )

  axis_list <- vector("list", num_axes)
  names(axis_list) <- paste0("axis", seq_len(num_axes))
  for (axis_idx in seq_len(num_axes)) {
    axis_list[[axis_idx]] <- matrix(result[, axis_idx + 2], nrow = nrows, ncol = ncols, byrow = TRUE)
  }

  append(output, axis_list)
}

popmaps_cell_coords <- function(raster_surface, legacy_compat = FALSE) {
  if (isTRUE(legacy_compat)) {
    return(popmaps_legacy_cell_coords(
      nrows = raster_surface@nrows,
      ncols = raster_surface@ncols,
      xmin = raster_surface@extent@xmin,
      ymax = raster_surface@extent@ymax,
      cell_size = raster::res(raster_surface)[1]
    ))
  }

  raster::xyFromCell(raster_surface, seq_len(raster::ncell(raster_surface)))
}

popmaps_legacy_cell_coords <- function(nrows, ncols, xmin, ymax, cell_size) {
  x_coords <- xmin + (cell_size / 2) + (seq_len(ncols) * cell_size)
  x_coords[ncols] <- xmin + (cell_size / 2) + ((ncols - 1) * cell_size)

  y_coords <- ymax + (cell_size / 2) - (seq_len(nrows) * cell_size)
  y_coords[nrows] <- ymax + (cell_size / 2) - ((nrows - 1) * cell_size)

  as.matrix(expand.grid(x = x_coords, y = y_coords))
}

popmaps_cell_site_distances <- function(coords, sampling_loc_coords) {
  distances <- matrix(NA_real_, nrow = nrow(coords), ncol = nrow(sampling_loc_coords))

  for (site_idx in seq_len(nrow(sampling_loc_coords))) {
    distances[, site_idx] <- popmaps_earth_dist(
      lat1 = coords[, 2],
      long1 = coords[, 1],
      lat2 = sampling_loc_coords[site_idx, 2],
      long2 = sampling_loc_coords[site_idx, 1]
    )
  }

  distances
}

popmaps_empirical_site_distances <- function(sampling_loc_coords) {
  num_emp_sites <- nrow(sampling_loc_coords)
  distances <- matrix(0, nrow = num_emp_sites, ncol = num_emp_sites)

  for (site_idx in seq_len(num_emp_sites)) {
    distances[, site_idx] <- popmaps_earth_dist(
      lat1 = sampling_loc_coords[, 2],
      long1 = sampling_loc_coords[, 1],
      lat2 = sampling_loc_coords[site_idx, 2],
      long2 = sampling_loc_coords[site_idx, 1]
    )
  }

  distances
}

popmaps_select_empirical_sites <- function(candidate_sites,
                                          empirical_distances,
                                          empirical_pt_dist,
                                          num_tested) {
  selected_sites <- integer(0)

  for (candidate_site in candidate_sites) {
    if (
      length(selected_sites) == 0 ||
        all(empirical_distances[candidate_site, selected_sites] >= empirical_pt_dist)
    ) {
      selected_sites <- c(selected_sites, candidate_site)
    }

    if (length(selected_sites) == num_tested) {
      return(selected_sites)
    }
  }

  stop(
    "`num_sites` did not include enough empirical sites separated by `empirical_pt_dist`.",
    call. = FALSE
  )
}

popmaps_earth_dist <- function(lat1, long1, lat2, long2) {
  rad <- pi / 180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlat <- b1 - a1
  dlon <- b2 - a2
  a <- (sin(dlat / 2))^2 + cos(a1) * cos(b1) * (sin(dlon / 2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145

  R * c
}
