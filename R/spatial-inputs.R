popmaps_prepare_inputs <- function(input_raster,
                                   input_locs,
                                   surface = c("G", "C"),
                                   num_sites = NULL,
                                   num_tested = NULL,
                                   threshold = NULL,
                                   ncore = NULL,
                                   empirical_pt_dist = NULL,
                                   popmod = NULL,
                                   jackknife = FALSE,
                                   require_legacy_c = TRUE) {
  surface <- match.arg(surface)

  if (surface == "C" && isTRUE(require_legacy_c) && !requireNamespace("gdistance", quietly = TRUE)) {
    stop("The 'gdistance' package is required when surface = 'C'.", call. = FALSE)
  }

  raster_input <- popmaps_prepare_raster(input_raster)
  locations <- popmaps_prepare_locations(input_locs)

  popmaps_validate_model_args(
    locations = locations,
    num_sites = num_sites,
    num_tested = num_tested,
    threshold = threshold,
    ncore = ncore,
    empirical_pt_dist = empirical_pt_dist,
    popmod = popmod,
    jackknife = jackknife
  )

  list(
    raster = raster_input$raster,
    rast = raster_input$rast,
    locations = locations,
    surface = surface
  )
}

popmaps_prepare_raster <- function(input_raster) {
  if (missing(input_raster) || identical(input_raster, "")) {
    stop("`input_raster` must be a raster object or a path to a raster file.", call. = FALSE)
  }

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
      stop("`input_raster` must be a terra SpatRaster, raster RasterLayer, or readable raster path.", call. = FALSE)
    }
  )

  if (terra::nlyr(rast) != 1) {
    stop("`input_raster` must contain exactly one raster layer.", call. = FALSE)
  }

  if (any(!is.finite(terra::res(rast))) || any(terra::res(rast) <= 0)) {
    stop("`input_raster` must have a finite, positive cell resolution.", call. = FALSE)
  }

  list(
    rast = rast,
    raster = raster::raster(rast)
  )
}

popmaps_prepare_locations <- function(input_locs) {
  if (missing(input_locs) || identical(input_locs, "")) {
    stop("`input_locs` must be a data frame or matrix.", call. = FALSE)
  }

  if (!is.data.frame(input_locs) && !is.matrix(input_locs)) {
    stop("`input_locs` must be a data frame or matrix.", call. = FALSE)
  }

  locations <- as.data.frame(input_locs, stringsAsFactors = FALSE)

  if (ncol(locations) < 4) {
    stop("`input_locs` must have at least four columns: site, longitude, latitude, and one ancestry coefficient.", call. = FALSE)
  }

  names(locations) <- paste0("V", seq_len(ncol(locations)))

  locations$V2 <- popmaps_numeric_column(locations$V2, "`input_locs` column 2 (longitude)")
  locations$V3 <- popmaps_numeric_column(locations$V3, "`input_locs` column 3 (latitude)")

  ancestry_cols <- seq.int(4, ncol(locations))
  for (idx in ancestry_cols) {
    locations[[idx]] <- popmaps_numeric_column(
      locations[[idx]],
      paste0("`input_locs` column ", idx, " (ancestry coefficient)")
    )
  }

  if (any(!is.finite(locations$V2)) || any(!is.finite(locations$V3))) {
    stop("`input_locs` longitude and latitude columns must contain finite values.", call. = FALSE)
  }

  ancestry <- as.matrix(locations[, ancestry_cols, drop = FALSE])
  if (any(!is.finite(ancestry))) {
    stop("`input_locs` ancestry coefficient columns must contain finite values.", call. = FALSE)
  }

  if (any(ancestry < 0)) {
    stop("`input_locs` ancestry coefficients must be non-negative.", call. = FALSE)
  }

  locations
}

popmaps_numeric_column <- function(x, label) {
  numeric_x <- suppressWarnings(as.numeric(x))

  if (any(is.na(numeric_x) & !is.na(x))) {
    stop(label, " must be numeric.", call. = FALSE)
  }

  numeric_x
}

popmaps_validate_model_args <- function(locations,
                                        num_sites = NULL,
                                        num_tested = NULL,
                                        threshold = NULL,
                                        ncore = NULL,
                                        empirical_pt_dist = NULL,
                                        popmod = NULL,
                                        jackknife = FALSE) {
  n_empirical <- nrow(locations)
  max_sites <- if (isTRUE(jackknife)) n_empirical - 1 else n_empirical

  if (!is.null(num_sites)) {
    popmaps_check_whole_number(num_sites, "`num_sites`")
    if (num_sites > max_sites) {
      stop("`num_sites` cannot exceed the number of available empirical sites.", call. = FALSE)
    }
  }

  if (!is.null(num_tested)) {
    popmaps_check_whole_number(num_tested, "`num_tested`")
    if (!is.null(num_sites) && num_tested > num_sites) {
      stop("`num_tested` cannot exceed `num_sites`.", call. = FALSE)
    }
    if (num_tested > max_sites) {
      stop("`num_tested` cannot exceed the number of available empirical sites.", call. = FALSE)
    }
  }

  if (!is.null(ncore)) {
    popmaps_check_whole_number(ncore, "`ncore`")
  }

  if (!is.null(threshold)) {
    popmaps_check_finite_scalar(threshold, "`threshold`")
  }

  if (!is.null(empirical_pt_dist)) {
    popmaps_check_finite_scalar(empirical_pt_dist, "`empirical_pt_dist`")
    if (empirical_pt_dist < 0) {
      stop("`empirical_pt_dist` must be non-negative.", call. = FALSE)
    }
  }

  if (!is.null(popmod)) {
    popmaps_check_finite_scalar(popmod, "`popmod`")
  }

  invisible(TRUE)
}

popmaps_check_whole_number <- function(x, label) {
  popmaps_check_finite_scalar(x, label)
  if (x < 1 || x != floor(x)) {
    stop(label, " must be a positive whole number.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_check_finite_scalar <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop(label, " must be one finite numeric value.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_prepare_point <- function(dec_long, dec_lat) {
  popmaps_check_finite_scalar(dec_long, "`dec_long`")
  popmaps_check_finite_scalar(dec_lat, "`dec_lat`")
  cbind(dec_long, dec_lat)
}
