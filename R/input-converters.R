#' Create a POPMAPS surface from coordinates
#'
#' @description
#' `surface_from_points()` builds a simple prediction grid around sampling
#' coordinates, then returns a [prepare_popmaps_surface()] object. This is useful
#' when users have empirical ancestry locations but have not supplied a raster
#' yet. The resulting surface is usually geographic (`surface = "G"`), although
#' a constant conductance surface can also be created with `surface = "C"`.
#'
#' @param points A data frame, matrix, `sf` object, or `terra::SpatVector`
#'   containing point coordinates. Data frames may be raw `x`/`y` coordinates or
#'   POPMAPS-style location tables where columns 2 and 3 are longitude and
#'   latitude.
#' @param coordinate_cols Optional coordinate column names or integer positions
#'   for data frame or matrix inputs. If `NULL`, common coordinate names are
#'   detected first; otherwise POPMAPS-style columns 2 and 3 are used when the
#'   first column is non-numeric, and columns 1 and 2 are used as a final
#'   fallback.
#' @param resolution Raster cell size in coordinate units. If `NULL`, a
#'   resolution is chosen so the longest coordinate span has about
#'   `target_cells` cells.
#' @param buffer Extra distance added around the coordinate bounding box. If
#'   `NULL`, the larger of one cell or five percent of the longest coordinate
#'   span is used.
#' @param target_cells Approximate number of cells along the longest coordinate
#'   axis when `resolution = NULL`.
#' @param crs Coordinate reference system assigned to the created raster. The
#'   default `""` leaves the CRS unset; pass a CRS string when the generated
#'   grid should carry projection metadata.
#' @param values Numeric scalar or vector of raster values. A scalar fills the
#'   grid. A vector must have one value per raster cell.
#' @inheritParams prepare_popmaps_surface
#'
#' @return A `popmaps_surface` object.
#'
#' @export
surface_from_points <- function(points,
                                coordinate_cols = NULL,
                                resolution = NULL,
                                buffer = NULL,
                                target_cells = 100,
                                crs = "",
                                values = 1,
                                surface = c("G", "C"),
                                surface_values = c("suitability", "conductance", "resistance"),
                                rescale_conductance = FALSE,
                                resistance_epsilon = sqrt(.Machine$double.eps)) {
  surface <- match.arg(surface)
  surface_values <- if (surface == "C") match.arg(surface_values) else "suitability"
  coords <- popmaps_coordinates_from_input(points, coordinate_cols = coordinate_cols, crs = crs)
  template <- popmaps_template_from_coords(
    coords = coords,
    resolution = resolution,
    buffer = buffer,
    target_cells = target_cells,
    crs = crs
  )

  template <- popmaps_fill_raster_values(template, values)
  prepare_popmaps_surface(
    input_raster = template,
    surface = surface,
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
}

#' Convert sf points to a POPMAPS location table
#'
#' @description
#' `locs_from_sf()` converts `sf` point features into the data frame shape
#' expected by `popmaps()`, `tune_popmaps()`, and
#' [compare_popmaps_surfaces()]. The returned table has site, longitude,
#' latitude, and ancestry coefficient columns.
#'
#' @param x An `sf` object with point geometries.
#' @param ancestry_cols Ancestry coefficient columns, supplied as names or
#'   integer positions. If `NULL`, all numeric non-geometry columns except
#'   `site_col` are used.
#' @param site_col Optional site identifier column, supplied as a name or integer
#'   position. If `NULL`, a column named `site`, `sample`, `id`, `name`, or
#'   `population` is used when present; otherwise synthetic site names are
#'   created.
#' @param crs CRS to transform the coordinates to before extracting x/y values.
#'   The default is `"EPSG:4326"` for longitude/latitude output. Set `crs = NULL`
#'   to keep the input CRS.
#' @param coord_names Names for the returned coordinate columns.
#'
#' @return A data frame accepted by POPMAPS modeling functions.
#'
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   locs <- locs_from_sf(sf::st_as_sf(
#'     data.frame(site = "a", axis1 = 1, axis2 = 0, lon = -110, lat = 39),
#'     coords = c("lon", "lat"),
#'     crs = 4326
#'   ))
#' }
#'
#' @export
locs_from_sf <- function(x,
                         ancestry_cols = NULL,
                         site_col = NULL,
                         crs = "EPSG:4326",
                         coord_names = c("lon", "lat")) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("The 'sf' package is required for `locs_from_sf()`.", call. = FALSE)
  }
  if (!inherits(x, "sf")) {
    stop("`x` must be an sf object.", call. = FALSE)
  }
  if (!is.character(coord_names) || length(coord_names) != 2 || anyNA(coord_names) || any(!nzchar(coord_names))) {
    stop("`coord_names` must contain two non-empty names.", call. = FALSE)
  }

  geometry_types <- unique(as.character(sf::st_geometry_type(x)))
  if (!all(geometry_types %in% c("POINT"))) {
    stop("`x` must contain POINT geometries.", call. = FALSE)
  }
  if (!is.null(crs)) {
    x <- sf::st_transform(x, crs)
  }

  coords <- sf::st_coordinates(x)
  if (ncol(coords) < 2 || any(!is.finite(coords[, 1:2]))) {
    stop("`x` must have finite point coordinates.", call. = FALSE)
  }

  attrs <- sf::st_drop_geometry(x)
  site_idx <- popmaps_resolve_site_column(attrs, site_col)
  site_values <- popmaps_extract_site_values(attrs, site_idx)
  ancestry_idx <- popmaps_resolve_ancestry_columns(attrs, ancestry_cols, site_idx)

  out <- data.frame(
    site = site_values,
    x = coords[, 1],
    y = coords[, 2],
    stringsAsFactors = FALSE
  )
  names(out)[2:3] <- coord_names
  out <- cbind(out, attrs[, ancestry_idx, drop = FALSE])
  popmaps_prepare_locations(out)
  out
}

#' Create candidate surfaces from a raster stack
#'
#' @description
#' `surfaces_from_raster_stack()` converts each layer of a multi-layer raster
#' into a named `popmaps_surface` object. This is a convenience helper for
#' comparing several user-supplied candidate surfaces with
#' [compare_popmaps_surfaces()].
#'
#' @param input_raster A multi-layer `terra::SpatRaster`, legacy raster stack or
#'   brick, or raster file path.
#' @param names Optional candidate surface names. Defaults to raster layer names
#'   when available.
#' @param include_geographic Logical. If `TRUE`, add a geographic surface using
#'   the first layer's geometry.
#' @param geographic_name Name for the optional geographic surface.
#' @inheritParams prepare_popmaps_surface
#'
#' @return A named list of `popmaps_surface` objects.
#'
#' @export
surfaces_from_raster_stack <- function(input_raster,
                                       surface = "C",
                                       surface_values = "suitability",
                                       names = NULL,
                                       include_geographic = FALSE,
                                       geographic_name = "geographic",
                                       mask = NULL,
                                       barrier = NULL,
                                       rescale_conductance = FALSE,
                                       resistance_epsilon = sqrt(.Machine$double.eps)) {
  rast <- popmaps_prepare_multilayer_raster(input_raster)
  n_layers <- terra::nlyr(rast)
  surface <- popmaps_recycle_surface_arg(surface, n_layers, "`surface`")
  surface_values <- popmaps_recycle_surface_arg(surface_values, n_layers, "`surface_values`")
  rescale_conductance <- popmaps_recycle_surface_arg(rescale_conductance, n_layers, "`rescale_conductance`")
  resistance_epsilon <- popmaps_recycle_surface_arg(resistance_epsilon, n_layers, "`resistance_epsilon`")

  surface_names <- popmaps_surface_stack_names(rast, names)
  surfaces <- lapply(seq_len(n_layers), function(idx) {
    prepare_popmaps_surface(
      input_raster = rast[[idx]],
      surface = surface[[idx]],
      surface_values = surface_values[[idx]],
      mask = mask,
      barrier = barrier,
      rescale_conductance = rescale_conductance[[idx]],
      resistance_epsilon = resistance_epsilon[[idx]]
    )
  })
  names(surfaces) <- surface_names

  if (isTRUE(include_geographic)) {
    if (!is.character(geographic_name) || length(geographic_name) != 1 || is.na(geographic_name) || !nzchar(geographic_name)) {
      stop("`geographic_name` must be one non-empty character value.", call. = FALSE)
    }
    if (geographic_name %in% names(surfaces)) {
      stop("`geographic_name` must not duplicate a candidate surface name.", call. = FALSE)
    }
    surfaces <- c(
      stats::setNames(
        list(prepare_popmaps_surface(rast[[1]], surface = "G", mask = mask)),
        geographic_name
      ),
      surfaces
    )
  }

  surfaces
}

#' Convert EEMS output to a POPMAPS conductance surface
#'
#' @description
#' `surface_from_eems()` converts raster-like EEMS outputs or coordinate/value
#' tables to a `surface = "C"` POPMAPS surface. Tables must contain x/y
#' coordinates plus one numeric surface-value column. High migration values
#' should usually be treated as conductance; effective resistance-style outputs
#' can be passed with `surface_values = "resistance"`.
#'
#' @param input A raster-like object/path, data frame, or matrix. Tabular inputs
#'   are rasterized using `template` or a grid derived from their coordinates.
#' @param value_col Surface-value column for tabular inputs. If `NULL`, common
#'   names such as `conductance`, `migration`, `value`, `mean`, `m`, or `w` are
#'   detected.
#' @inheritParams surface_from_points
#' @inheritParams prepare_popmaps_surface
#' @param template Optional raster defining output geometry for tabular inputs.
#'
#' @return A `popmaps_surface` object.
#'
#' @export
surface_from_eems <- function(input,
                              value_col = NULL,
                              coordinate_cols = NULL,
                              template = NULL,
                              resolution = NULL,
                              buffer = NULL,
                              target_cells = 100,
                              crs = "",
                              surface_values = c("conductance", "resistance", "suitability"),
                              rescale_conductance = FALSE,
                              resistance_epsilon = sqrt(.Machine$double.eps)) {
  surface_values <- match.arg(surface_values)
  popmaps_surface_from_external_values(
    input = input,
    source = "EEMS",
    value_col = value_col,
    coordinate_cols = coordinate_cols,
    template = template,
    resolution = resolution,
    buffer = buffer,
    target_cells = target_cells,
    crs = crs,
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
}

#' Convert FEEMS output to a POPMAPS conductance surface
#'
#' @description
#' `surface_from_feems()` converts raster-like FEEMS outputs or coordinate/value
#' tables to a `surface = "C"` POPMAPS surface. For tabular FEEMS exports, edge
#' or node weights where larger values imply easier movement should usually be
#' treated as conductance.
#'
#' @inheritParams surface_from_eems
#'
#' @return A `popmaps_surface` object.
#'
#' @export
surface_from_feems <- function(input,
                               value_col = NULL,
                               coordinate_cols = NULL,
                               template = NULL,
                               resolution = NULL,
                               buffer = NULL,
                               target_cells = 100,
                               crs = "",
                               surface_values = c("conductance", "resistance", "suitability"),
                               rescale_conductance = FALSE,
                               resistance_epsilon = sqrt(.Machine$double.eps)) {
  surface_values <- match.arg(surface_values)
  popmaps_surface_from_external_values(
    input = input,
    source = "FEEMS",
    value_col = value_col,
    coordinate_cols = coordinate_cols,
    template = template,
    resolution = resolution,
    buffer = buffer,
    target_cells = target_cells,
    crs = crs,
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
}

popmaps_coordinates_from_input <- function(points, coordinate_cols = NULL, crs = "EPSG:4326") {
  if (inherits(points, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("The 'sf' package is required for sf coordinate inputs.", call. = FALSE)
    }
    if (!is.null(crs) && nzchar(crs)) {
      points <- sf::st_transform(points, crs)
    }
    coords <- sf::st_coordinates(points)
    return(popmaps_validate_coordinate_matrix(coords[, 1:2, drop = FALSE]))
  }

  if (inherits(points, "SpatVector")) {
    if (!all(terra::geomtype(points) == "points")) {
      stop("`points` must contain point geometries.", call. = FALSE)
    }
    if (!is.null(crs) && nzchar(crs) && nzchar(terra::crs(points, describe = FALSE))) {
      points <- terra::project(points, crs)
    }
    coords <- terra::crds(points)
    return(popmaps_validate_coordinate_matrix(coords[, 1:2, drop = FALSE]))
  }

  if (!is.data.frame(points) && !is.matrix(points)) {
    stop("`points` must be a data frame, matrix, sf object, or terra SpatVector.", call. = FALSE)
  }

  data <- as.data.frame(points, stringsAsFactors = FALSE)
  coord_idx <- popmaps_resolve_coordinate_columns(data, coordinate_cols)
  coords <- data[, coord_idx, drop = FALSE]
  popmaps_validate_coordinate_matrix(coords)
}

popmaps_resolve_coordinate_columns <- function(data, coordinate_cols = NULL) {
  if (!is.null(coordinate_cols)) {
    if (is.character(coordinate_cols)) {
      if (!all(coordinate_cols %in% names(data))) {
        stop("`coordinate_cols` names must exist in the input data.", call. = FALSE)
      }
      return(match(coordinate_cols, names(data)))
    }
    if (is.numeric(coordinate_cols) && length(coordinate_cols) == 2) {
      if (any(coordinate_cols < 1) || any(coordinate_cols > ncol(data))) {
        stop("`coordinate_cols` positions are outside the input data.", call. = FALSE)
      }
      return(as.integer(coordinate_cols))
    }
    stop("`coordinate_cols` must be two column names or positions.", call. = FALSE)
  }

  lower_names <- tolower(names(data))
  candidates <- list(
    c("lon", "lat"),
    c("longitude", "latitude"),
    c("x", "y"),
    c("easting", "northing")
  )
  for (candidate in candidates) {
    if (all(candidate %in% lower_names)) {
      return(match(candidate, lower_names))
    }
  }

  if (ncol(data) >= 3 && !is.numeric(data[[1]])) {
    return(c(2L, 3L))
  }
  if (ncol(data) >= 2) {
    return(c(1L, 2L))
  }

  stop("Could not infer coordinate columns.", call. = FALSE)
}

popmaps_validate_coordinate_matrix <- function(coords) {
  coords <- as.matrix(coords[, 1:2, drop = FALSE])
  storage.mode(coords) <- "double"
  if (nrow(coords) < 1 || ncol(coords) != 2 || any(!is.finite(coords))) {
    stop("Coordinate inputs must contain finite x/y values.", call. = FALSE)
  }
  colnames(coords) <- c("x", "y")
  coords
}

popmaps_template_from_coords <- function(coords,
                                         resolution = NULL,
                                         buffer = NULL,
                                         target_cells = 100,
                                         crs = "") {
  popmaps_check_positive_optional(resolution, "`resolution`")
  popmaps_check_positive_optional(buffer, "`buffer`", allow_zero = TRUE)
  popmaps_check_whole_count(target_cells, "`target_cells`", positive = TRUE)

  x_range <- range(coords[, 1])
  y_range <- range(coords[, 2])
  span <- c(diff(x_range), diff(y_range))
  max_span <- max(span)

  if (is.null(resolution)) {
    if (max_span <= 0) {
      stop("`resolution` is required when all coordinates are identical.", call. = FALSE)
    }
    resolution <- max_span / target_cells
  }
  if (is.null(buffer)) {
    buffer <- max(resolution, max_span * 0.05)
  }

  xmin <- x_range[[1]] - buffer
  xmax <- x_range[[2]] + buffer
  ymin <- y_range[[1]] - buffer
  ymax <- y_range[[2]] + buffer

  ncols <- max(1L, ceiling((xmax - xmin) / resolution))
  nrows <- max(1L, ceiling((ymax - ymin) / resolution))
  xmax <- xmin + ncols * resolution
  ymax <- ymin + nrows * resolution

  terra::rast(
    nrows = nrows,
    ncols = ncols,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    crs = crs
  )
}

popmaps_fill_raster_values <- function(rast, values) {
  if (!is.numeric(values) || anyNA(values)) {
    stop("`values` must be numeric and cannot contain `NA`.", call. = FALSE)
  }
  n_cells <- terra::ncell(rast)
  if (length(values) == 1) {
    terra::values(rast) <- rep(values, n_cells)
  } else if (length(values) == n_cells) {
    terra::values(rast) <- values
  } else {
    stop("`values` must be one number or one value per raster cell.", call. = FALSE)
  }

  rast
}

popmaps_prepare_multilayer_raster <- function(input_raster) {
  rast <- tryCatch(
    {
      if (inherits(input_raster, "SpatRaster")) {
        input_raster
      } else if (inherits(input_raster, "Raster")) {
        terra::rast(input_raster)
      } else if (is.character(input_raster) && length(input_raster) == 1) {
        terra::rast(input_raster)
      } else {
        stop("unsupported", call. = FALSE)
      }
    },
    error = function(err) {
      stop("`input_raster` must be a terra SpatRaster, raster object, or readable raster path.", call. = FALSE)
    }
  )
  if (terra::nlyr(rast) < 1) {
    stop("`input_raster` must contain at least one raster layer.", call. = FALSE)
  }

  rast
}

popmaps_surface_stack_names <- function(rast, names = NULL) {
  n_layers <- terra::nlyr(rast)
  surface_names <- names
  if (is.null(surface_names)) {
    surface_names <- names(rast)
    if (is.null(surface_names) || any(is.na(surface_names)) || any(!nzchar(surface_names))) {
      surface_names <- paste0("surface_", seq_len(n_layers))
    }
  }
  if (!is.character(surface_names) || length(surface_names) != n_layers || anyNA(surface_names) || any(!nzchar(surface_names))) {
    stop("`names` must contain one non-empty name per raster layer.", call. = FALSE)
  }
  if (anyDuplicated(surface_names)) {
    stop("Candidate surface names must be unique.", call. = FALSE)
  }

  surface_names
}

popmaps_recycle_surface_arg <- function(x, n, label) {
  if (length(x) == 1) {
    return(rep(as.list(x), n))
  }
  if (length(x) != n) {
    stop(label, " must have length 1 or one value per raster layer.", call. = FALSE)
  }

  as.list(x)
}

popmaps_surface_from_external_values <- function(input,
                                                 source,
                                                 value_col,
                                                 coordinate_cols,
                                                 template,
                                                 resolution,
                                                 buffer,
                                                 target_cells,
                                                 crs,
                                                 surface_values,
                                                 rescale_conductance,
                                                 resistance_epsilon) {
  raster_surface <- popmaps_try_external_raster(input)
  if (!is.null(raster_surface)) {
    return(prepare_popmaps_surface(
      input_raster = raster_surface,
      surface = "C",
      surface_values = surface_values,
      rescale_conductance = rescale_conductance,
      resistance_epsilon = resistance_epsilon
    ))
  }

  values <- popmaps_read_external_surface_table(input, source = source)
  coord_idx <- popmaps_resolve_coordinate_columns(values, coordinate_cols)
  value_idx <- popmaps_resolve_value_column(values, value_col)
  coords <- popmaps_validate_coordinate_matrix(values[, coord_idx, drop = FALSE])
  surface_value <- suppressWarnings(as.numeric(values[[value_idx]]))
  if (any(is.na(surface_value) & !is.na(values[[value_idx]])) || any(!is.finite(surface_value))) {
    stop("Surface value column must contain finite numeric values.", call. = FALSE)
  }

  if (is.null(template)) {
    rast <- popmaps_template_from_coords(
      coords = coords,
      resolution = resolution,
      buffer = buffer,
      target_cells = target_cells,
      crs = crs
    )
  } else {
    rast <- popmaps_prepare_raster(template)$rast
  }

  point_values <- data.frame(x = coords[, 1], y = coords[, 2], value = surface_value)
  points <- terra::vect(point_values, geom = c("x", "y"), crs = crs)
  if (!is.na(terra::crs(rast, describe = FALSE)) &&
      nzchar(terra::crs(rast, describe = FALSE)) &&
      !is.na(terra::crs(points, describe = FALSE)) &&
      nzchar(terra::crs(points, describe = FALSE))) {
    points <- terra::project(points, terra::crs(rast))
  }
  surface_rast <- terra::rasterize(points, rast, field = "value", fun = mean, background = NA_real_)
  if (!any(!is.na(terra::values(surface_rast, mat = FALSE)))) {
    stop("No tabular surface values fell inside the output raster grid.", call. = FALSE)
  }

  prepare_popmaps_surface(
    input_raster = surface_rast,
    surface = "C",
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
}

popmaps_try_external_raster <- function(input) {
  if (inherits(input, "SpatRaster") || inherits(input, "Raster")) {
    return(input)
  }
  if (!is.character(input) || length(input) != 1 || !file.exists(input)) {
    return(NULL)
  }

  extension <- tolower(tools::file_ext(input))
  if (extension %in% c("csv", "tsv", "txt", "dat")) {
    return(NULL)
  }
  tryCatch(terra::rast(input), error = function(err) NULL)
}

popmaps_read_external_surface_table <- function(input, source) {
  if (is.data.frame(input) || is.matrix(input)) {
    return(as.data.frame(input, stringsAsFactors = FALSE))
  }
  if (!is.character(input) || length(input) != 1 || !file.exists(input)) {
    stop(source, " input must be a raster, data frame, matrix, or readable file path.", call. = FALSE)
  }

  extension <- tolower(tools::file_ext(input))
  if (extension == "csv") {
    return(utils::read.csv(input, stringsAsFactors = FALSE))
  }
  if (extension %in% c("tsv", "txt", "dat")) {
    return(utils::read.table(input, header = TRUE, stringsAsFactors = FALSE))
  }

  stop(source, " input path was not readable as a raster or supported table.", call. = FALSE)
}

popmaps_resolve_value_column <- function(data, value_col = NULL) {
  if (!is.null(value_col)) {
    if (is.character(value_col)) {
      if (!value_col %in% names(data)) {
        stop("`value_col` must name a column in the input data.", call. = FALSE)
      }
      return(match(value_col, names(data)))
    }
    if (is.numeric(value_col) && length(value_col) == 1 && value_col >= 1 && value_col <= ncol(data)) {
      return(as.integer(value_col))
    }
    stop("`value_col` must be one column name or position.", call. = FALSE)
  }

  lower_names <- tolower(names(data))
  candidate_names <- c(
    "conductance", "suitability", "resistance", "value", "rate",
    "migration", "edge_weight", "weight", "w", "m", "q", "mean", "posterior"
  )
  match_idx <- match(candidate_names, lower_names, nomatch = 0)
  match_idx <- match_idx[match_idx > 0]
  if (length(match_idx) > 0) {
    return(match_idx[[1]])
  }

  numeric_cols <- which(vapply(data, function(x) {
    numeric_x <- suppressWarnings(as.numeric(x))
    !any(is.na(numeric_x) & !is.na(x))
  }, logical(1)))
  if (length(numeric_cols) >= 3) {
    return(numeric_cols[[3]])
  }

  stop("Could not infer the surface value column; supply `value_col`.", call. = FALSE)
}

popmaps_resolve_site_column <- function(attrs, site_col) {
  if (!is.null(site_col)) {
    return(popmaps_resolve_single_column(attrs, site_col, "`site_col`"))
  }
  lower_names <- tolower(names(attrs))
  candidates <- c("site", "sample", "id", "name", "population")
  match_idx <- match(candidates, lower_names, nomatch = 0)
  match_idx <- match_idx[match_idx > 0]
  if (length(match_idx) > 0) {
    return(match_idx[[1]])
  }

  integer()
}

popmaps_extract_site_values <- function(attrs, site_idx) {
  if (length(site_idx) == 1) {
    return(as.character(attrs[[site_idx]]))
  }
  paste0("site_", seq_len(nrow(attrs)))
}

popmaps_resolve_ancestry_columns <- function(attrs, ancestry_cols, site_idx) {
  if (!is.null(ancestry_cols)) {
    idx <- popmaps_resolve_columns(attrs, ancestry_cols, "`ancestry_cols`")
  } else {
    idx <- setdiff(seq_along(attrs), site_idx)
    idx <- idx[vapply(attrs[idx], is.numeric, logical(1))]
  }
  if (length(idx) < 1) {
    stop("No ancestry coefficient columns were found.", call. = FALSE)
  }

  idx
}

popmaps_resolve_columns <- function(data, columns, label) {
  if (is.character(columns)) {
    if (!all(columns %in% names(data))) {
      stop(label, " names must exist in the input data.", call. = FALSE)
    }
    return(match(columns, names(data)))
  }
  if (is.numeric(columns)) {
    if (any(columns < 1) || any(columns > ncol(data))) {
      stop(label, " positions are outside the input data.", call. = FALSE)
    }
    return(as.integer(columns))
  }

  stop(label, " must be column names or positions.", call. = FALSE)
}

popmaps_resolve_single_column <- function(data, column, label) {
  idx <- popmaps_resolve_columns(data, column, label)
  if (length(idx) != 1) {
    stop(label, " must identify exactly one column.", call. = FALSE)
  }

  idx
}

popmaps_check_positive_optional <- function(x, label, allow_zero = FALSE) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }
  popmaps_check_finite_scalar(x, label)
  if (x < if (allow_zero) 0 else .Machine$double.eps) {
    stop(label, " must be positive.", call. = FALSE)
  }

  invisible(TRUE)
}
