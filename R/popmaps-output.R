#' Convert POPMAPS output to a terra raster
#'
#' @description Converts the list returned by [popmaps()] into a named
#' `terra::SpatRaster` using the extent, resolution, and coordinate reference
#' system from the original input raster.
#'
#' @param pop_raster_list A list returned by [popmaps()].
#' @param input_raster The raster used to create `pop_raster_list`. It may be a
#'   `terra::SpatRaster`, a legacy `raster::RasterLayer`, or a path readable by
#'   `terra::rast()`.
#' @param layer_names Optional character vector of layer names. By default,
#'   layers are named `hard_boundary`, `ancestry_probability`, `axis_1`,
#'   `axis_2`, and so on.
#'
#' @return A `terra::SpatRaster` with one layer per `popmaps()` output matrix.
#'
#' @examples
#' ex_raster <- raster::aggregate(hija_raster, fact = 240)
#' aps <- popmaps(
#'   input_raster = ex_raster,
#'   input_locs = hija_struc,
#'   surface = "G",
#'   empirical_pt_dist = 0,
#'   num_sites = 5,
#'   num_tested = 2,
#'   popmod = -0.05,
#'   threshold = 0,
#'   ncore = 1
#' )
#' popmaps_rast(aps, ex_raster)
#'
#' @export
popmaps_rast <- function(pop_raster_list, input_raster, layer_names = NULL) {
  if (!is.list(pop_raster_list) || length(pop_raster_list) < 3) {
    stop("`pop_raster_list` must be a list returned by popmaps().", call. = FALSE)
  }

  raster_input <- popmaps_prepare_raster(input_raster)
  template <- raster_input$rast
  expected_dim <- as.integer(c(terra::nrow(template), terra::ncol(template)))

  layer_names <- popmaps_output_layer_names(pop_raster_list, layer_names)
  values <- matrix(NA_real_, nrow = terra::ncell(template), ncol = length(pop_raster_list))

  for (idx in seq_along(pop_raster_list)) {
    surface <- popmaps_validate_output_matrix(
      pop_raster_list[[idx]],
      expected_dim = expected_dim,
      label = paste0("`pop_raster_list[[", idx, "]]`")
    )
    values[, idx] <- as.vector(t(surface))
  }

  output <- terra::rast(
    nrows = expected_dim[1],
    ncols = expected_dim[2],
    nlyrs = length(pop_raster_list),
    ext = terra::ext(template),
    crs = terra::crs(template)
  )
  names(output) <- layer_names
  terra::values(output) <- values

  output
}

#' Write POPMAPS output layers as GeoTIFFs
#'
#' @description Writes each layer returned by [popmaps_rast()] to a separate
#' GeoTIFF. The function returns a small manifest containing the layer names and
#' written file paths.
#'
#' @param pop_raster_list A list returned by [popmaps()].
#' @param input_raster The raster used to create `pop_raster_list`.
#' @param dir Output directory for GeoTIFF files. Created if it does not exist.
#' @param prefix File name prefix used for each layer.
#' @param layer_names Optional character vector of layer names passed to
#'   [popmaps_rast()].
#' @param overwrite Logical. If `FALSE`, the function stops before writing when
#'   any output file already exists.
#' @param ... Additional arguments passed to `terra::writeRaster()`.
#'
#' @return A data frame with `layer` and `path` columns.
#'
#' @examples
#' \dontrun{
#' paths <- write_popmaps(aps, ex_raster, dir = "outputs", prefix = "hija")
#' }
#'
#' @export
write_popmaps <- function(pop_raster_list,
                          input_raster,
                          dir = ".",
                          prefix = "popmaps",
                          layer_names = NULL,
                          overwrite = FALSE,
                          ...) {
  popmaps_check_character_scalar(dir, "`dir`")
  popmaps_check_character_scalar(prefix, "`prefix`")
  if (!nzchar(prefix)) {
    stop("`prefix` must not be empty.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }

  raster_output <- popmaps_rast(
    pop_raster_list = pop_raster_list,
    input_raster = input_raster,
    layer_names = layer_names
  )

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(dir)) {
    stop("Could not create output directory `dir`.", call. = FALSE)
  }

  file_names <- paste0(prefix, "_", names(raster_output), ".tif")
  paths <- file.path(dir, file_names)
  existing_paths <- paths[file.exists(paths)]
  if (!overwrite && length(existing_paths)) {
    stop(
      "Output file(s) already exist. Use `overwrite = TRUE` to replace them: ",
      paste(existing_paths, collapse = ", "),
      call. = FALSE
    )
  }

  for (idx in seq_along(paths)) {
    terra::writeRaster(raster_output[[idx]], paths[idx], overwrite = overwrite, ...)
  }

  data.frame(
    layer = names(raster_output),
    path = normalizePath(paths, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

popmaps_output_layer_names <- function(pop_raster_list, layer_names = NULL) {
  default_names <- c(
    "hard_boundary",
    "ancestry_probability",
    paste0("axis_", seq_len(length(pop_raster_list) - 2))
  )

  if (is.null(layer_names)) {
    return(default_names)
  }

  if (!is.character(layer_names) ||
      length(layer_names) != length(pop_raster_list) ||
      any(is.na(layer_names)) ||
      any(!nzchar(layer_names))) {
    stop(
      "`layer_names` must be a non-empty character vector with one name per output layer.",
      call. = FALSE
    )
  }

  make.names(layer_names, unique = TRUE)
}

popmaps_validate_output_matrix <- function(surface, expected_dim, label) {
  if (!is.matrix(surface) && !is.data.frame(surface)) {
    stop(label, " must be a matrix-like object.", call. = FALSE)
  }

  surface <- as.matrix(surface)
  if (!all(dim(surface) == expected_dim)) {
    stop(
      label,
      " dimensions must match `input_raster` rows and columns.",
      call. = FALSE
    )
  }

  storage.mode(surface) <- "double"
  surface
}

popmaps_check_character_scalar <- function(x, label) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop(label, " must be one character value.", call. = FALSE)
  }

  invisible(TRUE)
}
