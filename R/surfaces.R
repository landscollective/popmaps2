#' Prepare a POPMAPS interpolation surface
#'
#' @description
#' `prepare_popmaps_surface()` records the scientific meaning of a raster before
#' it is used for POPMAPS interpolation or surface comparison. Geographic
#' surfaces (`surface = "G"`) use only the raster geometry. Suitability- or
#' conductance-weighted surfaces (`surface = "C"`) use raster values to define
#' relative movement or gene-flow connectivity.
#'
#' This helper creates a validated surface object for modern least-cost
#' interpolation, tuning, and candidate-surface comparison workflows.
#'
#' @param input_raster A `terra::SpatRaster`, legacy `raster::RasterLayer`, or
#'   path to a raster file.
#' @param surface Character. `"G"` uses raster geometry only. `"C"` uses raster
#'   values as a conductance-like landscape surface.
#' @param surface_values Character. Meaning of raster values when
#'   `surface = "C"`. `"suitability"` and `"conductance"` use values directly;
#'   `"resistance"` converts values to conductance using an inverse transform.
#'   Ignored when `surface = "G"`.
#' @param mask Optional raster with the same geometry as `input_raster`.
#'   Non-zero, non-`NA` cells mark cells eligible for prediction. This is kept
#'   separate from conductance because a prediction mask is not necessarily a
#'   movement barrier.
#' @param barrier Optional raster with the same geometry as `input_raster`.
#'   Non-zero, non-`NA` cells mark cells that should be treated as
#'   non-traversable in modern least-cost workflows.
#' @param rescale_conductance Logical. If `TRUE`, rescale non-missing
#'   conductance values to the range 0-1 after any resistance conversion. The
#'   default is `FALSE` to preserve legacy POPMAPS behavior for MaxEnt logistic
#'   suitability rasters.
#' @param resistance_epsilon Positive numeric scalar added to resistance values
#'   before inversion to avoid infinite conductance when resistance is zero.
#'
#' @return A `popmaps_surface` object containing the original grid, optional
#'   conductance raster, optional mask and barrier rasters, and transformation
#'   metadata.
#'
#' @examples
#' geographic <- prepare_popmaps_surface(hija_raster, surface = "G")
#' suitability <- prepare_popmaps_surface(
#'   hija_raster,
#'   surface = "C",
#'   surface_values = "suitability"
#' )
#'
#' @export
prepare_popmaps_surface <- function(input_raster,
                                    surface = c("G", "C"),
                                    surface_values = c("suitability", "conductance", "resistance"),
                                    mask = NULL,
                                    barrier = NULL,
                                    rescale_conductance = FALSE,
                                    resistance_epsilon = sqrt(.Machine$double.eps)) {
  surface <- match.arg(surface)
  surface_values_supplied <- !missing(surface_values)

  if (surface == "C") {
    surface_values <- match.arg(surface_values)
  } else {
    if (surface_values_supplied) {
      surface_values <- match.arg(surface_values)
    }
    surface_values <- NA_character_
  }

  if (!is.logical(rescale_conductance) || length(rescale_conductance) != 1 || is.na(rescale_conductance)) {
    stop("`rescale_conductance` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  raster_input <- popmaps_prepare_raster(input_raster)
  rast <- raster_input$rast
  names(rast) <- "surface"

  mask_rast <- popmaps_prepare_optional_surface_raster(mask, rast, "`mask`")
  barrier_rast <- popmaps_prepare_optional_surface_raster(barrier, rast, "`barrier`")

  conductance <- NULL
  transform <- "geometry_only"

  if (surface == "C") {
    conductance <- popmaps_prepare_conductance(
      rast = rast,
      surface_values = surface_values,
      rescale_conductance = rescale_conductance,
      resistance_epsilon = resistance_epsilon
    )

    if (!is.null(barrier_rast)) {
      conductance_values <- terra::values(conductance, mat = FALSE)
      barrier_values <- terra::values(barrier_rast, mat = FALSE)
      conductance_values[barrier_values == 1] <- NA_real_
      terra::values(conductance) <- conductance_values
      popmaps_assert_positive_conductance(conductance)
    }

    transform <- attr(conductance, "popmaps_transform", exact = TRUE)
    attr(conductance, "popmaps_transform") <- NULL
  }

  structure(
    list(
      rast = rast,
      raster = raster_input$raster,
      surface = surface,
      surface_values = surface_values,
      conductance = conductance,
      mask = mask_rast,
      barrier = barrier_rast,
      rescale_conductance = rescale_conductance,
      transform = transform,
      resistance_epsilon = if (identical(surface_values, "resistance")) resistance_epsilon else NA_real_
    ),
    class = "popmaps_surface"
  )
}

#' @export
print.popmaps_surface <- function(x, ...) {
  cat("<popmaps_surface>\n")
  cat("  surface: ", x$surface, "\n", sep = "")

  if (identical(x$surface, "C")) {
    cat("  surface values: ", x$surface_values, "\n", sep = "")
    cat("  transform: ", x$transform, "\n", sep = "")
    cat("  rescaled: ", x$rescale_conductance, "\n", sep = "")
  } else {
    cat("  raster values: ignored for distances\n")
  }

  cat("  mask: ", if (is.null(x$mask)) "none" else "provided", "\n", sep = "")
  cat("  barrier: ", if (is.null(x$barrier)) "none" else "provided", "\n", sep = "")
  invisible(x)
}

popmaps_prepare_conductance <- function(rast,
                                        surface_values,
                                        rescale_conductance,
                                        resistance_epsilon) {
  raster_values <- terra::values(rast, mat = FALSE)
  observed <- !is.na(raster_values)

  if (any(!is.finite(raster_values[observed]))) {
    stop("`input_raster` values must be finite or `NA`.", call. = FALSE)
  }

  if (any(raster_values[observed] < 0)) {
    stop("`input_raster` values must be non-negative for `surface = \"C\"`.", call. = FALSE)
  }

  conductance <- rast
  transform <- "identity"

  if (surface_values == "resistance") {
    popmaps_check_finite_scalar(resistance_epsilon, "`resistance_epsilon`")
    if (resistance_epsilon <= 0) {
      stop("`resistance_epsilon` must be positive.", call. = FALSE)
    }

    conductance_values <- rep(NA_real_, length(raster_values))
    conductance_values[observed] <- 1 / (raster_values[observed] + resistance_epsilon)
    terra::values(conductance) <- conductance_values
    transform <- "inverse_resistance"
  }

  popmaps_assert_positive_conductance(conductance)

  if (isTRUE(rescale_conductance)) {
    conductance <- popmaps_rescale_conductance(conductance)
    transform <- paste(transform, "rescaled_0_1", sep = "+")
  }

  names(conductance) <- "conductance"
  attr(conductance, "popmaps_transform") <- transform
  conductance
}

popmaps_assert_positive_conductance <- function(conductance) {
  conductance_values <- terra::values(conductance, mat = FALSE)
  conductance_observed <- !is.na(conductance_values)
  if (!any(conductance_values[conductance_observed] > 0)) {
    stop("`surface = \"C\"` requires at least one positive conductance value.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_rescale_conductance <- function(conductance) {
  conductance_values <- terra::values(conductance, mat = FALSE)
  observed <- !is.na(conductance_values)
  conductance_min <- min(conductance_values[observed])
  conductance_max <- max(conductance_values[observed])

  if (conductance_max == conductance_min) {
    conductance_values[observed] <- 1
  } else {
    conductance_values[observed] <- (conductance_values[observed] - conductance_min) /
      (conductance_max - conductance_min)
  }

  terra::values(conductance) <- conductance_values
  conductance
}

popmaps_prepare_optional_surface_raster <- function(input_raster, template, label) {
  if (is.null(input_raster)) {
    return(NULL)
  }

  prepared <- popmaps_prepare_raster(input_raster)$rast
  if (!isTRUE(terra::compareGeom(prepared, template, stopOnError = FALSE))) {
    stop(label, " must have the same geometry as `input_raster`.", call. = FALSE)
  }

  values <- terra::values(prepared, mat = FALSE)
  observed <- !is.na(values)
  if (any(!is.finite(values[observed]))) {
    stop(label, " values must be finite or `NA`.", call. = FALSE)
  }

  bool_values <- as.numeric(observed & values != 0)
  terra::values(prepared) <- bool_values
  names(prepared) <- gsub("`", "", label)
  prepared
}
