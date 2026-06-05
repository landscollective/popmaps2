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
                                 resistance_epsilon,
                                 legacy_compat = FALSE,
                                 quiet = TRUE) {
  nrows <- raster_surface@nrows
  ncols <- raster_surface@ncols

  num_axes <- ncol(species_data) - 3
  ancestry <- as.matrix(species_data[, seq.int(4, ncol(species_data)), drop = FALSE])
  sampling_loc_coords <- as.matrix(species_data[, 2:3, drop = FALSE])

  coords <- popmaps_cell_coords(raster_surface, legacy_compat = legacy_compat)
  raster_values <- raster::extract(raster_surface, coords)
  geographic_cell_distances <- popmaps_cell_site_distances(coords, sampling_loc_coords)

  popmaps_inform("Preparing conductance surface and graph.", quiet = quiet)
  surface_object <- prepare_popmaps_surface(
    input_raster = raster_surface,
    surface = "C",
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
  graph <- popmaps_cost_distance_graph(surface_object, directions = 8)

  popmaps_inform("Calculating least-cost distances among empirical sites.", quiet = quiet)
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
  popmaps_inform("Calculating least-cost distances from empirical sites to raster cells.", quiet = quiet)
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
    popmaps_progress_message(cell_idx, nrow(coords), "Estimated ancestry for raster cell", quiet = quiet)
    if (is.na(raster_values[cell_idx])) {
      next
    }

    if (isTRUE(legacy_compat)) {
      geographic_order <- order(geographic_cell_distances[cell_idx, ])[seq_len(num_sites)]
      cost_values <- cost_cell_distances[cell_idx, geographic_order]
      cell_order <- geographic_order[order(cost_values, na.last = TRUE)]
    } else {
      reachable_sites <- which(is.finite(cost_cell_distances[cell_idx, ]))
      if (length(reachable_sites) < num_sites) {
        next
      }
      cell_order <- reachable_sites[order(cost_cell_distances[cell_idx, reachable_sites])][seq_len(num_sites)]
    }
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
