#' Validate against the POPMAPS 1.03 baseline
#'
#' @description Runs a small Hilaria jamesii ancestry probability surface and
#' compares the result to a frozen POPMAPS 1.03 reference output bundled with
#' the package. This is intended as a quick scientific regression check before
#' deeper optimization work changes the modeling internals.
#'
#' @param tolerance Numeric tolerance used when comparing current output to the
#'   frozen reference output.
#' @param quiet Logical. If `FALSE`, prints a short pass/fail message.
#'
#' @return A data frame with one row per output surface and columns describing
#'   maximum absolute difference, mean absolute difference, NA mismatches, and
#'   pass/fail status. The overall result is stored in the `passed` attribute.
#'
#' @examples
#' validate_popmaps_baseline()
#'
#' @export
validate_popmaps_baseline <- function(tolerance = sqrt(.Machine$double.eps),
                                      quiet = FALSE) {
  popmaps_check_finite_scalar(tolerance, "`tolerance`")
  if (tolerance < 0) {
    stop("`tolerance` must be non-negative.", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1 || is.na(quiet)) {
    stop("`quiet` must be TRUE or FALSE.", call. = FALSE)
  }

  reference_path <- system.file(
    "extdata",
    "validation",
    "popmaps-1.03-hija-small-reference.rds",
    package = "popmaps2",
    mustWork = TRUE
  )
  reference <- readRDS(reference_path)
  params <- reference$params

  validation_data <- new.env(parent = emptyenv())
  utils::data(
    list = c("hija_raster", "hija_struc"),
    package = "popmaps2",
    envir = validation_data
  )

  ex_raster <- raster::aggregate(validation_data$hija_raster, fact = params$aggregate_fact)

  current <- popmaps(
    input_raster = ex_raster,
    input_locs = validation_data$hija_struc,
    surface = params$surface,
    empirical_pt_dist = params$empirical_pt_dist,
    num_sites = params$num_sites,
    num_tested = params$num_tested,
    popmod = params$popmod,
    threshold = params$threshold,
    ncore = 1
  )

  comparison <- popmaps_compare_output(
    current = current,
    reference = reference$result,
    tolerance = tolerance
  )
  passed <- all(comparison$passed)
  attr(comparison, "passed") <- passed
  attr(comparison, "params") <- params
  attr(comparison, "reference") <- reference$metadata

  if (!quiet) {
    if (passed) {
      message("POPMAPS 1.03 baseline validation passed.")
    } else {
      message("POPMAPS 1.03 baseline validation failed.")
    }
  }

  comparison
}

popmaps_compare_output <- function(current, reference, tolerance) {
  if (!is.list(current) || !is.list(reference) || length(current) != length(reference)) {
    stop("Current and reference outputs must be lists with the same length.", call. = FALSE)
  }

  surface_names <- c(
    "hard_boundary",
    "ancestry_probability",
    paste0("axis_", seq_len(length(current) - 2))
  )

  rows <- vector("list", length(current))
  for (idx in seq_along(current)) {
    current_surface <- current[[idx]]
    reference_surface <- reference[[idx]]
    same_dim <- identical(dim(current_surface), dim(reference_surface))

    if (same_dim) {
      current_values <- as.numeric(current_surface)
      reference_values <- as.numeric(reference_surface)
      na_mismatch <- sum(is.na(current_values) != is.na(reference_values))
      diffs <- abs(current_values - reference_values)
      finite_diffs <- diffs[is.finite(diffs)]
      max_abs_diff <- if (length(finite_diffs)) max(finite_diffs) else 0
      mean_abs_diff <- if (length(finite_diffs)) mean(finite_diffs) else 0
      passed <- na_mismatch == 0 && max_abs_diff <= tolerance
      n_cell <- length(current_values)
    } else {
      na_mismatch <- NA_integer_
      max_abs_diff <- Inf
      mean_abs_diff <- Inf
      passed <- FALSE
      n_cell <- NA_integer_
    }

    rows[[idx]] <- data.frame(
      surface = surface_names[idx],
      n_cell = n_cell,
      max_abs_diff = max_abs_diff,
      mean_abs_diff = mean_abs_diff,
      na_mismatch = na_mismatch,
      passed = passed,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}
