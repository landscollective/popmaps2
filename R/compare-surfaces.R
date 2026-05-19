#' Compare candidate POPMAPS interpolation surfaces
#'
#' @description
#' `compare_popmaps_surfaces()` runs matched `tune_popmaps()` validation across
#' multiple candidate surfaces. It asks which user-supplied surface best predicts
#' withheld empirical ancestry estimates in the POPMAPS interpolation workflow.
#' This is predictive model selection for ancestry surfaces, not a replacement
#' for upstream SDM, resistance-surface, EEMS/FEEMS, or landscape-genetic
#' hypothesis testing.
#'
#' @param input_locs A data frame or matrix with sampling location name,
#'   longitude, latitude, and one or more ancestry coefficient columns.
#' @param surfaces Named list of candidate surfaces. Each element may be a
#'   `popmaps_surface` object returned by [prepare_popmaps_surface()] or a list
#'   with `input_raster` or `raster`, `surface`, and optional
#'   `surface_values`, `mask`, `barrier`, `rescale_conductance`, and
#'   `resistance_epsilon` entries.
#' @param near_best_tolerance Non-negative relative tolerance for labeling
#'   surfaces as statistically near-best. The default keeps surfaces within 5%
#'   of the best primary metric score.
#' @inheritParams tune_popmaps
#'
#' @return A `popmaps_surface_comparison` object with:
#' \describe{
#'   \item{summary}{One row per surface, ranked by the primary metric.}
#'   \item{best}{The best-ranked surface row.}
#'   \item{near_best}{Surfaces within `near_best_tolerance` of the best score.}
#'   \item{support}{One-row conservative interpretation of surface support.}
#'   \item{tunings}{Named list of `popmaps_tuning` objects, one per surface.}
#' }
#'
#' @examples
#' ex_raster <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
#' terra::values(ex_raster) <- 1
#' locs <- data.frame(
#'   site = paste0("s", 1:4),
#'   lon = c(0.5, 3.5, 0.5, 3.5),
#'   lat = c(0.5, 0.5, 3.5, 3.5),
#'   axis1 = c(0.9, 0.8, 0.2, 0.1),
#'   axis2 = c(0.1, 0.2, 0.8, 0.9)
#' )
#' surfaces <- list(
#'   geographic = prepare_popmaps_surface(ex_raster, surface = "G"),
#'   suitability = prepare_popmaps_surface(
#'     ex_raster,
#'     surface = "C",
#'     surface_values = "suitability"
#'   )
#' )
#' comparison <- compare_popmaps_surfaces(
#'   input_locs = locs,
#'   surfaces = surfaces,
#'   empirical_pt_dist = 0,
#'   num_sites = 3,
#'   num_tested = 2,
#'   popmod = -0.1,
#'   quiet = TRUE
#' )
#' comparison$summary
#'
#' @export
compare_popmaps_surfaces <- function(input_locs,
                                     surfaces,
                                     empirical_pt_dist = 5,
                                     num_sites = 10,
                                     num_tested = c(2, 3, 4, 5, 6, 7, 8),
                                     popmod = c(-0.001, -0.01, -0.05, -0.1, -0.15),
                                     threshold = 0,
                                     validation = c("loo", "spatial_block"),
                                     n_blocks = 4,
                                     block_assignments = NULL,
                                     spatial_block_repeats = 1,
                                     spatial_block_seed = NULL,
                                     primary_metric = c("rmse", "mae", "hellinger",
                                                        "dominant_accuracy",
                                                        "dominant_probability"),
                                     dist_prob_func = function(popmod_temp, distance) {
                                       exp(popmod_temp * distance)
                                     },
                                     near_best_tolerance = 0.05,
                                     quiet = TRUE) {
  validation <- match.arg(validation)
  primary_metric <- match.arg(primary_metric)
  if (!is.function(dist_prob_func)) {
    stop("`dist_prob_func` must be a function.", call. = FALSE)
  }
  popmaps_check_finite_scalar(near_best_tolerance, "`near_best_tolerance`")
  if (near_best_tolerance < 0) {
    stop("`near_best_tolerance` must be non-negative.", call. = FALSE)
  }

  surface_specs <- popmaps_normalize_surface_specs(surfaces)
  if (
    validation == "spatial_block" &&
      is.null(block_assignments) &&
      spatial_block_repeats > 1 &&
      is.null(spatial_block_seed)
  ) {
    spatial_block_seed <- 1L
  }

  tunings <- vector("list", length(surface_specs))
  names(tunings) <- names(surface_specs)
  for (surface_name in names(surface_specs)) {
    surface_object <- surface_specs[[surface_name]]
    tunings[[surface_name]] <- tune_popmaps(
      input_raster = surface_object,
      input_locs = input_locs,
      surface = surface_object$surface,
      empirical_pt_dist = empirical_pt_dist,
      num_sites = num_sites,
      num_tested = num_tested,
      popmod = popmod,
      threshold = threshold,
      validation = validation,
      n_blocks = n_blocks,
      block_assignments = block_assignments,
      spatial_block_repeats = spatial_block_repeats,
      spatial_block_seed = spatial_block_seed,
      primary_metric = primary_metric,
      dist_prob_func = dist_prob_func,
      quiet = TRUE
    )
  }

  summary <- popmaps_surface_comparison_summary(
    tunings = tunings,
    surfaces = surface_specs,
    primary_metric = primary_metric
  )
  best <- popmaps_surface_comparison_best(summary, primary_metric)
  near_best <- popmaps_surface_comparison_near_best(
    summary = summary,
    primary_metric = primary_metric,
    near_best_tolerance = near_best_tolerance
  )
  support <- popmaps_surface_comparison_support(
    summary = summary,
    best = best,
    near_best = near_best,
    primary_metric = primary_metric,
    near_best_tolerance = near_best_tolerance
  )

  comparison <- list(
    summary = summary,
    best = best,
    near_best = near_best,
    support = support,
    tunings = tunings,
    primary_metric = primary_metric,
    validation = validation,
    near_best_tolerance = near_best_tolerance,
    spatial_block_seed = spatial_block_seed,
    call = match.call()
  )
  class(comparison) <- "popmaps_surface_comparison"

  if (!isTRUE(quiet)) {
    message(
      "Compared ",
      nrow(summary),
      " candidate surfaces with ",
      validation,
      " validation."
    )
  }

  comparison
}

#' @export
print.popmaps_surface_comparison <- function(x, ...) {
  cat("POPMAPS surface comparison\n")
  cat("Primary metric: ", x$primary_metric, "\n", sep = "")
  cat("Validation: ", x$validation, "\n", sep = "")
  cat("Decision: ", x$support$decision, "\n\n", sep = "")
  cat("Surface ranking:\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

popmaps_normalize_surface_specs <- function(surfaces) {
  if (!is.list(surfaces) || length(surfaces) < 2) {
    stop("`surfaces` must be a named list with at least two candidate surfaces.", call. = FALSE)
  }

  surface_names <- names(surfaces)
  if (is.null(surface_names) || any(!nzchar(surface_names))) {
    surface_names <- paste0("surface_", seq_along(surfaces))
  }
  if (anyDuplicated(surface_names)) {
    stop("`surfaces` names must be unique.", call. = FALSE)
  }

  normalized <- lapply(seq_along(surfaces), function(idx) {
    popmaps_normalize_surface_spec(surfaces[[idx]], surface_names[[idx]])
  })
  names(normalized) <- surface_names
  normalized
}

popmaps_normalize_surface_spec <- function(spec, surface_name) {
  if (inherits(spec, "popmaps_surface")) {
    return(spec)
  }
  if (!is.list(spec)) {
    stop(
      "`surfaces[['", surface_name, "']]` must be a `popmaps_surface` or a surface specification list.",
      call. = FALSE
    )
  }

  input_raster <- if ("input_raster" %in% names(spec)) spec$input_raster else spec$raster
  if (is.null(input_raster)) {
    stop(
      "`surfaces[['", surface_name, "']]` must include `input_raster` or `raster`.",
      call. = FALSE
    )
  }

  surface <- if ("surface" %in% names(spec)) spec$surface else "G"
  surface_values <- if ("surface_values" %in% names(spec)) spec$surface_values else "suitability"
  mask <- if ("mask" %in% names(spec)) spec$mask else NULL
  barrier <- if ("barrier" %in% names(spec)) spec$barrier else NULL
  rescale_conductance <- if ("rescale_conductance" %in% names(spec)) spec$rescale_conductance else FALSE
  resistance_epsilon <- if ("resistance_epsilon" %in% names(spec)) {
    spec$resistance_epsilon
  } else {
    sqrt(.Machine$double.eps)
  }

  prepare_popmaps_surface(
    input_raster = input_raster,
    surface = surface,
    surface_values = surface_values,
    mask = mask,
    barrier = barrier,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
}

popmaps_surface_comparison_summary <- function(tunings, surfaces, primary_metric) {
  rows <- lapply(names(tunings), function(surface_name) {
    tuning <- tunings[[surface_name]]
    surface <- surfaces[[surface_name]]
    best <- tuning$best

    data.frame(
      surface_name = surface_name,
      surface = surface$surface,
      surface_values = surface$surface_values,
      transform = surface$transform,
      rescale_conductance = surface$rescale_conductance,
      primary_metric = primary_metric,
      score = best[[primary_metric]][1],
      distance_units = best$distance_units[1],
      num_sites = best$num_sites[1],
      num_tested = best$num_tested[1],
      popmod = best$popmod[1],
      half_distance = best$half_distance_km[1],
      ten_pct_distance = best$ten_pct_distance_km[1],
      empirical_pt_dist = best$empirical_pt_dist[1],
      n_combinations = nrow(tuning$results),
      n_scored = best$n_scored[1],
      failed_folds = best$failed_folds[1],
      n_validation_repeats = best$n_validation_repeats[1],
      n_validation_folds = best$n_validation_folds[1],
      stringsAsFactors = FALSE
    )
  })

  summary <- do.call(rbind, rows)
  summary <- summary[popmaps_order_surface_summary(summary, primary_metric), , drop = FALSE]
  rownames(summary) <- NULL
  summary$rank <- seq_len(nrow(summary))

  best_score <- summary$score[1]
  if (popmaps_metric_is_maximized(primary_metric)) {
    summary$delta_from_best <- best_score - summary$score
  } else {
    summary$delta_from_best <- summary$score - best_score
  }
  summary$percent_from_best <- vapply(
    summary$delta_from_best,
    popmaps_percent_change,
    numeric(1),
    reference = best_score
  )
  summary
}

popmaps_order_surface_summary <- function(summary, primary_metric) {
  maximize <- popmaps_metric_is_maximized(primary_metric)
  score <- summary$score
  score_order <- if (maximize) -score else score
  score_order[!is.finite(score_order)] <- Inf

  order(
    summary$failed_folds,
    -summary$n_scored,
    score_order,
    summary$surface_name
  )
}

popmaps_surface_comparison_best <- function(summary, primary_metric) {
  summary[popmaps_order_surface_summary(summary, primary_metric)[1], , drop = FALSE]
}

popmaps_surface_comparison_near_best <- function(summary,
                                                 primary_metric,
                                                 near_best_tolerance) {
  maximize <- popmaps_metric_is_maximized(primary_metric)
  best_score <- summary$score[1]
  metric_scale <- max(abs(best_score), .Machine$double.eps)
  scored <- is.finite(summary$score)
  near_best <- if (maximize) {
    scored & summary$score >= best_score - near_best_tolerance * metric_scale
  } else {
    scored & summary$score <= best_score + near_best_tolerance * metric_scale
  }

  summary[near_best, , drop = FALSE]
}

popmaps_surface_comparison_support <- function(summary,
                                               best,
                                               near_best,
                                               primary_metric,
                                               near_best_tolerance) {
  finite_scores <- is.finite(summary$score)
  decision <- if (!any(finite_scores)) {
    "insufficient_validation"
  } else if (best$failed_folds[1] > 0) {
    "unstable_or_incomplete"
  } else if (nrow(near_best) > 1) {
    "surfaces_indistinguishable"
  } else {
    paste0(best$surface_name[1], "_supported")
  }

  data.frame(
    decision = decision,
    best_surface = best$surface_name[1],
    primary_metric = primary_metric,
    best_score = best$score[1],
    near_best_tolerance = near_best_tolerance,
    n_near_best = nrow(near_best),
    near_best_surfaces = paste(near_best$surface_name, collapse = ", "),
    stringsAsFactors = FALSE
  )
}
