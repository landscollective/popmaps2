popmaps_cost_surface <- function(raster_surface,
                                 species_data,
                                 empirical_pt_dist,
                                 num_sites,
                                 num_tested,
                                 popmod,
                                 threshold,
                                 dist_prob_func,
                                 surface_values,
                                 rescale_conductance,
                                 resistance_epsilon) {
  nrows <- raster_surface@nrows
  ncols <- raster_surface@ncols
  ymax <- raster_surface@extent@ymax
  xmin <- raster_surface@extent@xmin
  cell_size <- raster::res(raster_surface)[1]

  num_axes <- ncol(species_data) - 3
  ancestry <- as.matrix(species_data[, seq.int(4, ncol(species_data)), drop = FALSE])
  sampling_loc_coords <- as.matrix(species_data[, 2:3, drop = FALSE])

  coords <- popmaps_legacy_cell_coords(
    nrows = nrows,
    ncols = ncols,
    xmin = xmin,
    ymax = ymax,
    cell_size = cell_size
  )
  raster_values <- raster::extract(raster_surface, coords)
  geographic_cell_distances <- popmaps_cell_site_distances(coords, sampling_loc_coords)

  surface_object <- prepare_popmaps_surface(
    input_raster = raster_surface,
    surface = "C",
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
  graph <- popmaps_cost_distance_graph(surface_object, directions = 8)

  empirical_distances <- tryCatch(
    popmaps_cost_distance_matrix(
      surface = surface_object,
      from_coords = sampling_loc_coords,
      directions = 8,
      graph = graph
    ),
    error = function(err) {
      stop(
        "Could not calculate least-cost distances among empirical sites: ",
        conditionMessage(err),
        call. = FALSE
      )
    }
  )
  cell_distance_lookup <- tryCatch(
    popmaps_cost_distance_to_cells(
      surface = surface_object,
      from_coords = sampling_loc_coords,
      directions = 8,
      graph = graph
    ),
    error = function(err) {
      stop(
        "Could not calculate least-cost distances from empirical sites to raster cells: ",
        conditionMessage(err),
        call. = FALSE
      )
    }
  )

  legacy_cells <- raster::cellFromXY(graph$raster, coords)
  cost_cell_distances <- matrix(
    NA_real_,
    nrow = nrow(coords),
    ncol = nrow(sampling_loc_coords)
  )
  valid_cells <- which(!is.na(raster_values) & !is.na(legacy_cells))
  cost_cell_distances[valid_cells, ] <- cell_distance_lookup[legacy_cells[valid_cells], , drop = FALSE]

  hard_boundaries <- max.col(ancestry, ties.method = "first")
  result <- matrix(NA_real_, nrow = nrow(coords), ncol = num_axes + 2)

  for (cell_idx in seq_len(nrow(coords))) {
    if (is.na(raster_values[cell_idx])) {
      next
    }

    geographic_order <- order(geographic_cell_distances[cell_idx, ])[seq_len(num_sites)]
    cost_values <- cost_cell_distances[cell_idx, geographic_order]
    cell_order <- geographic_order[order(cost_values, na.last = TRUE)]
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
      cost_cell_distances[cell_idx, selected_sites],
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

  axis_list <- vector("list", num_axes)
  names(axis_list) <- paste0("axis", seq_len(num_axes))
  for (axis_idx in seq_len(num_axes)) {
    axis_list[[axis_idx]] <- matrix(result[, axis_idx + 2], nrow = nrows, ncol = ncols, byrow = TRUE)
  }

  append(output, axis_list)
}
