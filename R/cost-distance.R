popmaps_cost_distance_graph <- function(surface, directions = 8) {
  if (!inherits(surface, "popmaps_surface")) {
    stop("`surface` must be a `popmaps_surface` object.", call. = FALSE)
  }
  if (!identical(surface$surface, "C") || is.null(surface$conductance)) {
    stop("Least-cost distances require a `surface = \"C\"` popmaps surface.", call. = FALSE)
  }
  if (!directions %in% c(4, 8)) {
    stop("`directions` must be 4 or 8.", call. = FALSE)
  }

  conductance <- surface$conductance
  conductance_values <- terra::values(conductance, mat = FALSE)
  traversable <- which(!is.na(conductance_values))
  positive <- traversable[conductance_values[traversable] > 0]

  if (length(positive) == 0) {
    stop("Least-cost distances require at least one positive conductance cell.", call. = FALSE)
  }

  raster_template <- raster::raster(conductance)
  adjacency <- popmaps_cost_adjacency(
    raster_template = raster_template,
    conductance_values = conductance_values,
    directions = directions
  )

  vertices <- data.frame(name = popmaps_cell_names(traversable), stringsAsFactors = FALSE)
  if (nrow(adjacency) == 0) {
    graph <- igraph::make_empty_graph(n = nrow(vertices), directed = FALSE)
    graph <- igraph::set_vertex_attr(graph, "name", value = vertices$name)
  } else {
    graph <- igraph::graph_from_data_frame(
      d = data.frame(
        from = popmaps_cell_names(adjacency$from),
        to = popmaps_cell_names(adjacency$to),
        weight = adjacency$cost,
        stringsAsFactors = FALSE
      ),
      directed = FALSE,
      vertices = vertices
    )
  }

  list(
    graph = graph,
    cells = traversable,
    raster = raster_template,
    directions = directions
  )
}

popmaps_cost_distance_matrix <- function(surface,
                                         from_coords,
                                         to_coords = NULL,
                                         directions = 8,
                                         graph = NULL) {
  if (is.null(to_coords)) {
    to_coords <- from_coords
  }

  if (is.null(graph)) {
    graph <- popmaps_cost_distance_graph(surface, directions = directions)
  }

  from_coords <- popmaps_coords_to_matrix(from_coords, "`from_coords`")
  to_coords <- popmaps_coords_to_matrix(to_coords, "`to_coords`")
  from_cells <- popmaps_coords_to_traversable_cells(graph, from_coords, "`from_coords`")
  to_cells <- popmaps_coords_to_traversable_cells(graph, to_coords, "`to_coords`")

  popmaps_igraph_distances(
    graph = graph$graph,
    from_cells = from_cells,
    to_cells = to_cells
  )
}

popmaps_cost_distance_to_cells <- function(surface,
                                           from_coords,
                                           directions = 8,
                                           graph = NULL) {
  if (is.null(graph)) {
    graph <- popmaps_cost_distance_graph(surface, directions = directions)
  }

  from_coords <- popmaps_coords_to_matrix(from_coords, "`from_coords`")
  from_cells <- popmaps_coords_to_traversable_cells(graph, from_coords, "`from_coords`")
  cell_distances <- matrix(
    NA_real_,
    nrow = raster::ncell(graph$raster),
    ncol = length(from_cells)
  )

  distances <- popmaps_igraph_distances(
    graph = graph$graph,
    from_cells = from_cells,
    to_cells = graph$cells
  )

  cell_distances[graph$cells, ] <- t(distances)
  colnames(cell_distances) <- rownames(from_coords)
  cell_distances
}

popmaps_new_distance_cache <- function(enabled = TRUE) {
  cache <- new.env(parent = emptyenv())
  cache$.enabled <- isTRUE(enabled)
  cache$.hits <- 0L
  cache$.misses <- 0L
  cache$.sets <- 0L
  cache
}

popmaps_distance_cache_get <- function(cache, key) {
  if (is.null(cache) || !isTRUE(cache$.enabled)) {
    return(NULL)
  }

  if (exists(key, envir = cache, inherits = FALSE)) {
    cache$.hits <- cache$.hits + 1L
    return(get(key, envir = cache, inherits = FALSE))
  }

  cache$.misses <- cache$.misses + 1L
  NULL
}

popmaps_distance_cache_set <- function(cache, key, value) {
  if (!is.null(cache) && isTRUE(cache$.enabled)) {
    assign(key, value, envir = cache)
    cache$.sets <- cache$.sets + 1L
  }

  value
}

popmaps_distance_cache_summary <- function(cache) {
  if (is.null(cache)) {
    return(data.frame(
      enabled = FALSE,
      entries = 0L,
      hits = 0L,
      misses = 0L,
      sets = 0L
    ))
  }

  keys <- setdiff(ls(envir = cache, all.names = TRUE), c(".enabled", ".hits", ".misses", ".sets"))
  data.frame(
    enabled = isTRUE(cache$.enabled),
    entries = length(keys),
    hits = cache$.hits,
    misses = cache$.misses,
    sets = cache$.sets
  )
}

popmaps_distance_cache_key <- function(...) {
  parts <- vapply(list(...), popmaps_distance_cache_part, character(1))
  paste(parts, collapse = "||")
}

popmaps_distance_cache_part <- function(x) {
  x <- paste(as.character(x), collapse = ",")
  gsub("[|[:space:]]+", "_", x)
}

popmaps_coord_signature <- function(coords) {
  coords <- popmaps_coords_to_matrix(coords, "`coords`")
  values <- c(
    nrow(coords),
    ncol(coords),
    colSums(coords),
    colSums(coords^2),
    apply(coords, 2, min),
    apply(coords, 2, max)
  )

  paste(signif(values, 12), collapse = ",")
}

popmaps_cached_empirical_site_distances <- function(coords,
                                                    distance_cache = NULL,
                                                    cache_key = "geographic") {
  key <- popmaps_distance_cache_key(
    "site_distances",
    "G",
    cache_key,
    popmaps_coord_signature(coords)
  )
  cached <- popmaps_distance_cache_get(distance_cache, key)
  if (!is.null(cached)) {
    return(cached)
  }

  popmaps_distance_cache_set(
    distance_cache,
    key,
    popmaps_empirical_site_distances(coords)
  )
}

popmaps_cached_cost_graph <- function(surface,
                                      directions = 8,
                                      distance_cache = NULL,
                                      cache_key = "surface") {
  key <- popmaps_distance_cache_key("cost_graph", cache_key, directions)
  cached <- popmaps_distance_cache_get(distance_cache, key)
  if (!is.null(cached)) {
    return(cached)
  }

  popmaps_distance_cache_set(
    distance_cache,
    key,
    popmaps_cost_distance_graph(surface, directions = directions)
  )
}

popmaps_cached_cost_site_distances <- function(surface,
                                               coords,
                                               directions = 8,
                                               graph = NULL,
                                               distance_cache = NULL,
                                               cache_key = "surface") {
  key <- popmaps_distance_cache_key(
    "site_distances",
    "C",
    cache_key,
    directions,
    popmaps_coord_signature(coords)
  )
  cached <- popmaps_distance_cache_get(distance_cache, key)
  if (!is.null(cached)) {
    return(cached)
  }

  if (is.null(graph)) {
    graph <- popmaps_cached_cost_graph(
      surface = surface,
      directions = directions,
      distance_cache = distance_cache,
      cache_key = cache_key
    )
  }

  popmaps_distance_cache_set(
    distance_cache,
    key,
    popmaps_cost_distance_matrix(
      surface = surface,
      from_coords = coords,
      directions = directions,
      graph = graph
    )
  )
}

popmaps_cost_adjacency <- function(raster_template,
                                   conductance_values,
                                   directions) {
  ncells <- raster::ncell(raster_template)
  nrows <- raster::nrow(raster_template)
  ncols <- raster::ncol(raster_template)
  cells <- seq_len(ncells)
  rows <- ((cells - 1) %/% ncols) + 1
  cols <- ((cells - 1) %% ncols) + 1

  offsets <- if (directions == 4) {
    data.frame(row = c(0, 1), col = c(1, 0))
  } else {
    data.frame(row = c(0, 1, 1, 1), col = c(1, -1, 0, 1))
  }

  from <- integer(0)
  to <- integer(0)
  for (idx in seq_len(nrow(offsets))) {
    to_rows <- rows + offsets$row[idx]
    to_cols <- cols + offsets$col[idx]
    valid <- to_rows >= 1 & to_rows <= nrows & to_cols >= 1 & to_cols <= ncols

    from <- c(from, cells[valid])
    to <- c(to, ((to_rows[valid] - 1) * ncols) + to_cols[valid])
  }

  edge_conductance <- rowMeans(cbind(conductance_values[from], conductance_values[to]))
  keep <- is.finite(edge_conductance) & edge_conductance > 0
  from <- from[keep]
  to <- to[keep]
  edge_conductance <- edge_conductance[keep]

  if (length(from) == 0) {
    return(data.frame(from = integer(0), to = integer(0), cost = numeric(0)))
  }

  edge_distance <- popmaps_edge_distances(
    raster_template = raster_template,
    from = from,
    to = to
  )

  data.frame(
    from = from,
    to = to,
    cost = edge_distance / edge_conductance
  )
}

popmaps_edge_distances <- function(raster_template, from, to) {
  from_xy <- raster::xyFromCell(raster_template, from)
  to_xy <- raster::xyFromCell(raster_template, to)
  edge_distance <- raster::pointDistance(
    from_xy,
    to_xy,
    longlat = raster::isLonLat(raster_template)
  )
  scale_value <- popmaps_geo_correction_scale(raster_template)

  as.numeric(edge_distance) / scale_value
}

popmaps_geo_correction_scale <- function(raster_template) {
  midpoint <- c(
    mean(c(raster::xmin(raster_template), raster::xmax(raster_template))),
    mean(c(raster::ymin(raster_template), raster::ymax(raster_template)))
  )
  scale_value <- raster::pointDistance(
    midpoint,
    midpoint + c(raster::xres(raster_template), 0),
    longlat = raster::isLonLat(raster_template)
  )
  as.numeric(scale_value)
}

popmaps_coords_to_matrix <- function(coords, label) {
  if (!is.matrix(coords) && is.numeric(coords) && length(coords) == 2) {
    coords <- matrix(coords, nrow = 1)
  }
  if (inherits(coords, "SpatialPoints")) {
    coords <- sp::coordinates(coords)
  }
  if (!is.matrix(coords) && !is.data.frame(coords)) {
    stop(label, " must be a two-column coordinate matrix or data frame.", call. = FALSE)
  }

  coords <- as.matrix(coords)
  if (ncol(coords) != 2) {
    stop(label, " must have exactly two columns.", call. = FALSE)
  }
  storage.mode(coords) <- "double"
  if (any(!is.finite(coords))) {
    stop(label, " must contain finite coordinates.", call. = FALSE)
  }

  coords
}

popmaps_coords_to_traversable_cells <- function(graph, coords, label) {
  cells <- raster::cellFromXY(graph$raster, coords)
  if (any(is.na(cells))) {
    stop(label, " must fall inside the interpolation raster.", call. = FALSE)
  }

  graph_cells <- as.integer(igraph::V(graph$graph)$name)
  if (any(!cells %in% graph_cells)) {
    stop(label, " must fall on traversable conductance cells.", call. = FALSE)
  }

  cells
}

popmaps_igraph_distances <- function(graph, from_cells, to_cells) {
  distances <- igraph::distances(
    graph,
    v = popmaps_cell_names(unique(from_cells)),
    to = popmaps_cell_names(unique(to_cells)),
    weights = igraph::E(graph)$weight,
    algorithm = "dijkstra"
  )

  distances[
    match(popmaps_cell_names(from_cells), rownames(distances)),
    match(popmaps_cell_names(to_cells), colnames(distances)),
    drop = FALSE
  ]
}

popmaps_cell_names <- function(cells) {
  format(as.integer(cells), scientific = FALSE, trim = TRUE)
}
