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
#' @param surface_grid Character. `"surface_specific"` derives `popmod` and
#'   `empirical_pt_dist` from distances measured over each candidate surface
#'   when either argument is `NULL`. `"shared"` derives missing values from
#'   geographic sampling-site distances and uses them for every surface.
#' @param near_best_tolerance Non-negative relative tolerance for labeling
#'   surfaces as statistically near-best. The default keeps surfaces within 5%
#'   of the best primary metric score.
#' @param cache Logical. If `TRUE`, reuse in-memory geographic and least-cost
#'   distance objects while resolving tuning grids and evaluating candidate
#'   surfaces. This is useful when several surfaces or validation designs reuse
#'   the same empirical coordinates.
#' @inheritParams tune_popmaps
#'
#' @return A `popmaps_surface_comparison` object with:
#' \describe{
#'   \item{summary}{One row per surface, ranked by the primary metric.}
#'   \item{best}{The best-ranked surface row.}
#'   \item{near_best}{Surfaces within `near_best_tolerance` of the best score.}
#'   \item{support}{One-row conservative interpretation of surface support.}
#'   \item{tunings}{Named list of `popmaps_tuning` objects, one per surface.}
#'   \item{grids}{Named list of surface-specific tuning grids used for the comparison.}
#'   \item{cache}{Summary of in-memory distance-cache use during the comparison.}
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
                                     empirical_pt_dist = NULL,
                                     num_sites = NULL,
                                     num_tested = NULL,
                                     popmod = NULL,
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
                                     surface_grid = c("surface_specific", "shared"),
                                     near_best_tolerance = 0.05,
                                     cache = TRUE,
                                     quiet = TRUE) {
  validation <- match.arg(validation)
  primary_metric <- match.arg(primary_metric)
  surface_grid <- match.arg(surface_grid)
  if (!is.function(dist_prob_func)) {
    stop("`dist_prob_func` must be a function.", call. = FALSE)
  }
  popmaps_check_finite_scalar(near_best_tolerance, "`near_best_tolerance`")
  if (near_best_tolerance < 0) {
    stop("`near_best_tolerance` must be non-negative.", call. = FALSE)
  }
  popmaps_check_logical_scalar(cache, "`cache`")

  surface_specs <- popmaps_normalize_surface_specs(surfaces)
  distance_cache <- popmaps_new_distance_cache(enabled = cache)
  if (
    validation == "spatial_block" &&
      is.null(block_assignments) &&
      spatial_block_repeats > 1 &&
      is.null(spatial_block_seed)
  ) {
    spatial_block_seed <- 1L
  }

  tunings <- vector("list", length(surface_specs))
  grids <- vector("list", length(surface_specs))
  names(tunings) <- names(surface_specs)
  names(grids) <- names(surface_specs)
  for (surface_name in names(surface_specs)) {
    surface_object <- surface_specs[[surface_name]]
    grids[[surface_name]] <- popmaps_resolve_comparison_grid(
      input_locs = input_locs,
      surface_object = surface_object,
      num_sites = num_sites,
      num_tested = num_tested,
      popmod = popmod,
      empirical_pt_dist = empirical_pt_dist,
      surface_grid = surface_grid,
      distance_cache = distance_cache,
      surface_cache_key = surface_name
    )
    parameter_grid <- popmaps_make_tuning_grid(
      num_sites = grids[[surface_name]]$num_sites,
      num_tested = grids[[surface_name]]$num_tested,
      popmod = grids[[surface_name]]$popmod,
      empirical_pt_dist = grids[[surface_name]]$empirical_pt_dist
    )
    tunings[[surface_name]] <- popmaps_evaluate_tuning_grid(
      input_raster = surface_object,
      input_locs = input_locs,
      surface = surface_object$surface,
      parameter_grid = parameter_grid,
      threshold = threshold,
      validation = validation,
      n_blocks = n_blocks,
      block_assignments = block_assignments,
      spatial_block_repeats = spatial_block_repeats,
      spatial_block_seed = spatial_block_seed,
      primary_metric = primary_metric,
      dist_prob_func = dist_prob_func,
      surface_values = surface_object$surface_values,
      rescale_conductance = surface_object$rescale_conductance,
      resistance_epsilon = surface_object$resistance_epsilon,
      quiet = TRUE,
      call = match.call(),
      distance_cache = distance_cache,
      surface_cache_key = surface_name
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
    grids = grids,
    primary_metric = primary_metric,
    validation = validation,
    surface_grid = surface_grid,
    near_best_tolerance = near_best_tolerance,
    cache = popmaps_distance_cache_summary(distance_cache),
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

popmaps_resolve_comparison_grid <- function(input_locs,
                                            surface_object,
                                            num_sites,
                                            num_tested,
                                            popmod,
                                            empirical_pt_dist,
                                            surface_grid,
                                            distance_cache = NULL,
                                            surface_cache_key = NULL) {
  suggested <- if (surface_grid == "surface_specific") {
    popmaps_suggest_surface_tuning_grid_cached(
      input_locs = input_locs,
      num_sites = num_sites,
      num_tested = num_tested,
      surface_object = surface_object,
      distance_cache = distance_cache,
      surface_cache_key = surface_cache_key
    )
  } else {
    suggest_tuning_grid(
      input_locs = input_locs,
      num_sites = num_sites,
      num_tested = num_tested
    )
  }

  suggested$popmod <- if (is.null(popmod)) {
    suggested$popmod
  } else {
    popmaps_check_tuning_values(popmod, "`popmod`")
  }
  suggested$empirical_pt_dist <- if (is.null(empirical_pt_dist)) {
    suggested$empirical_pt_dist
  } else {
    popmaps_check_tuning_values(
      empirical_pt_dist,
      "`empirical_pt_dist`",
      nonnegative = TRUE
    )
  }

  suggested
}

popmaps_suggest_surface_tuning_grid_cached <- function(input_locs,
                                                       surface_object,
                                                       num_sites,
                                                       num_tested,
                                                       distance_cache = NULL,
                                                       surface_cache_key = NULL,
                                                       empirical_pt_dist_probs = c(0.05, 0.10, 0.25),
                                                       distance_weights = c(0.95, 0.75, 0.50, 0.25, 0.10, 0.05),
                                                       distance_reference = "median",
                                                       max_num_tested = 8) {
  locations <- popmaps_prepare_locations(input_locs)
  coords <- as.matrix(locations[, 2:3, drop = FALSE])
  if (is.null(surface_cache_key)) {
    surface_cache_key <- surface_object$surface
  }

  if (identical(surface_object$surface, "G")) {
    distance_matrix <- popmaps_cached_empirical_site_distances(
      coords = coords,
      distance_cache = distance_cache,
      cache_key = surface_cache_key
    )
    surface_values <- NA_character_
    distance_units <- "km"
  } else {
    graph <- popmaps_cached_cost_graph(
      surface = surface_object,
      directions = 8,
      distance_cache = distance_cache,
      cache_key = surface_cache_key
    )
    distance_matrix <- popmaps_cached_cost_site_distances(
      surface = surface_object,
      coords = coords,
      directions = 8,
      graph = graph,
      distance_cache = distance_cache,
      cache_key = surface_cache_key
    )
    surface_values <- surface_object$surface_values
    distance_units <- "cost_distance"
  }

  popmaps_suggest_tuning_grid_from_distances(
    locations = locations,
    distance_matrix = distance_matrix,
    num_sites = num_sites,
    num_tested = num_tested,
    empirical_pt_dist_probs = empirical_pt_dist_probs,
    distance_weights = distance_weights,
    distance_reference = distance_reference,
    max_num_tested = max_num_tested,
    surface = surface_object$surface,
    surface_values = surface_values,
    distance_units = distance_units
  )
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

#' Plot a POPMAPS surface comparison
#'
#' @description
#' `plot_surface_comparison()` draws compact base R diagnostic plots from a
#' [compare_popmaps_surfaces()] result. These plots are intended to make surface
#' support interpretable without requiring users to inspect the comparison list
#' by hand.
#'
#' @param comparison A `popmaps_surface_comparison` object returned by
#'   [compare_popmaps_surfaces()].
#' @param type Plot type. `"score"` plots the primary validation metric by
#'   surface. `"percent_from_best"` plots relative loss from the best surface.
#'   `"best_parameters"` plots the best parameter values selected for each
#'   surface. `"score_distribution"` shows whether each surface has a sharp or
#'   broad tuning optimum. `"near_best_parameters"` plots full and near-best
#'   parameter ranges for each surface.
#' @param col Optional vector of plotting colors. Named vectors are matched to
#'   surface names; unnamed vectors are recycled in comparison order.
#' @param main Optional plot title.
#' @param ... Additional arguments passed to base plotting functions.
#'
#' @return Invisibly returns `comparison`.
#'
#' @export
plot_surface_comparison <- function(comparison,
                                    type = c("score", "percent_from_best", "best_parameters",
                                             "score_distribution", "near_best_parameters"),
                                    col = NULL,
                                    main = NULL,
                                    ...) {
  popmaps_check_surface_comparison(comparison)
  type <- match.arg(type)

  summary <- comparison$summary
  surface_names <- summary$surface_name
  colors <- popmaps_surface_comparison_colors(surface_names, col)

  if (type == "score") {
    popmaps_plot_surface_scores(summary, comparison$primary_metric, colors, main, ...)
  } else if (type == "percent_from_best") {
    popmaps_plot_surface_percent_from_best(summary, colors, main, ...)
  } else if (type == "best_parameters") {
    popmaps_plot_surface_best_parameters(summary, colors, main, ...)
  } else if (type == "score_distribution") {
    popmaps_plot_surface_score_distribution(comparison, colors, main, ...)
  } else {
    popmaps_plot_surface_near_best_parameters(comparison, colors, main, ...)
  }

  invisible(comparison)
}

#' Write a POPMAPS surface comparison report
#'
#' @description
#' `write_surface_comparison_report()` writes the main tables, diagnostic plots,
#' and a small Markdown interpretation file from a
#' [compare_popmaps_surfaces()] result. It is designed for reproducible local
#' validation runs where users need durable artifacts to compare candidate
#' geographic, suitability, conductance, or resistance surfaces.
#'
#' @inheritParams plot_surface_comparison
#' @param dir Output directory. Created recursively if needed.
#' @param prefix File prefix used for all report artifacts.
#' @param include_tuning_results Logical. If `TRUE`, write one full
#'   `tune_popmaps()` results CSV per candidate surface.
#' @param overwrite Logical. If `FALSE`, stop before replacing an existing
#'   output file.
#' @param width,height,res PNG device settings for report figures.
#'
#' @return A `popmaps_surface_comparison_report` list containing paths to the
#'   written report, tables, figures, and optional tuning-result CSVs.
#'
#' @export
write_surface_comparison_report <- function(comparison,
                                            dir = ".",
                                            prefix = "surface-comparison",
                                            include_tuning_results = TRUE,
                                            overwrite = TRUE,
                                            width = 1500,
                                            height = 900,
                                            res = 150) {
  popmaps_check_surface_comparison(comparison)
  popmaps_check_character_scalar(dir, "`dir`")
  popmaps_check_character_scalar(prefix, "`prefix`")
  popmaps_check_logical_scalar(include_tuning_results, "`include_tuning_results`")
  popmaps_check_logical_scalar(overwrite, "`overwrite`")
  popmaps_check_png_dimension(width, "`width`")
  popmaps_check_png_dimension(height, "`height`")
  popmaps_check_png_dimension(res, "`res`")

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir, mustWork = TRUE)
  figure_dir <- file.path(dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  table_paths <- c(
    summary = file.path(dir, paste0(prefix, "-summary.csv")),
    support = file.path(dir, paste0(prefix, "-support.csv")),
    grids = file.path(dir, paste0(prefix, "-grids.csv")),
    tuning_diagnostics = file.path(dir, paste0(prefix, "-tuning-diagnostics.csv")),
    near_best_parameters = file.path(dir, paste0(prefix, "-near-best-parameters.csv"))
  )
  figure_paths <- c(
    score = file.path(figure_dir, paste0(prefix, "-scores.png")),
    percent_from_best = file.path(figure_dir, paste0(prefix, "-percent-from-best.png")),
    best_parameters = file.path(figure_dir, paste0(prefix, "-best-parameters.png")),
    score_distribution = file.path(figure_dir, paste0(prefix, "-score-distribution.png")),
    near_best_parameters = file.path(figure_dir, paste0(prefix, "-near-best-parameters.png"))
  )
  report_path <- file.path(dir, paste0(prefix, "-report.md"))

  popmaps_check_report_paths(
    c(table_paths, figure_paths, report = report_path),
    overwrite = overwrite
  )

  grids <- popmaps_flatten_surface_comparison_grids(comparison$grids)
  diagnostics <- popmaps_surface_tuning_diagnostics(comparison)
  near_best_parameters <- popmaps_surface_near_best_parameter_ranges(comparison)
  popmaps_write_report_table(comparison$summary, table_paths[["summary"]])
  popmaps_write_report_table(comparison$support, table_paths[["support"]])
  popmaps_write_report_table(grids, table_paths[["grids"]])
  popmaps_write_report_table(diagnostics, table_paths[["tuning_diagnostics"]])
  popmaps_write_report_table(near_best_parameters, table_paths[["near_best_parameters"]])

  popmaps_write_surface_comparison_png(
    figure_paths[["score"]],
    width = width,
    height = height,
    res = res,
    expr = plot_surface_comparison(comparison, type = "score")
  )
  popmaps_write_surface_comparison_png(
    figure_paths[["percent_from_best"]],
    width = width,
    height = height,
    res = res,
    expr = plot_surface_comparison(comparison, type = "percent_from_best")
  )
  popmaps_write_surface_comparison_png(
    figure_paths[["best_parameters"]],
    width = width,
    height = height,
    res = res,
    expr = plot_surface_comparison(comparison, type = "best_parameters")
  )
  popmaps_write_surface_comparison_png(
    figure_paths[["score_distribution"]],
    width = width,
    height = height,
    res = res,
    expr = plot_surface_comparison(comparison, type = "score_distribution")
  )
  popmaps_write_surface_comparison_png(
    figure_paths[["near_best_parameters"]],
    width = width,
    height = height,
    res = res,
    expr = plot_surface_comparison(comparison, type = "near_best_parameters")
  )

  tuning_paths <- character()
  if (isTRUE(include_tuning_results)) {
    tuning_paths <- popmaps_write_surface_tuning_results(
      comparison = comparison,
      dir = dir,
      prefix = prefix,
      overwrite = overwrite
    )
  }

  report <- popmaps_surface_comparison_markdown(
    comparison = comparison,
    table_paths = table_paths,
    figure_paths = figure_paths,
    tuning_paths = tuning_paths
  )
  writeLines(report, report_path)

  manifest <- list(
    report = report_path,
    tables = table_paths,
    figures = figure_paths,
    tuning_results = tuning_paths
  )
  class(manifest) <- "popmaps_surface_comparison_report"
  manifest
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
      half_distance = best$half_distance[1],
      ten_pct_distance = best$ten_pct_distance[1],
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

popmaps_check_surface_comparison <- function(comparison) {
  if (!inherits(comparison, "popmaps_surface_comparison") || is.null(comparison$summary)) {
    stop(
      "`comparison` must be an object returned by `compare_popmaps_surfaces()`.",
      call. = FALSE
    )
  }

  required <- c("surface_name", "score", "percent_from_best")
  missing <- setdiff(required, names(comparison$summary))
  if (length(missing) > 0) {
    stop("`comparison$summary` is missing required columns.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_check_logical_scalar <- function(x, label) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(label, " must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_check_png_dimension <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop(label, " must be one positive finite number.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_surface_comparison_colors <- function(surface_names, col = NULL) {
  if (is.null(col)) {
    colors <- grDevices::hcl.colors(length(surface_names), "Dark 3")
  } else {
    if (!is.character(col) || length(col) < 1 || anyNA(col)) {
      stop("`col` must be a character vector of colors.", call. = FALSE)
    }
    if (!is.null(names(col)) && all(surface_names %in% names(col))) {
      colors <- unname(col[surface_names])
    } else {
      colors <- rep(col, length.out = length(surface_names))
    }
  }

  stats::setNames(colors, surface_names)
}

popmaps_surface_metric_label <- function(primary_metric) {
  goal <- if (popmaps_metric_is_maximized(primary_metric)) {
    "higher is better"
  } else {
    "lower is better"
  }

  paste0(primary_metric, " (", goal, ")")
}

popmaps_surface_plot_labels <- function(surface_names, width = 13) {
  vapply(surface_names, function(surface_name) {
    label <- gsub("[_]+", " ", surface_name)
    paste(strwrap(label, width = width), collapse = "\n")
  }, character(1))
}

popmaps_surface_tuning_diagnostics <- function(comparison) {
  popmaps_check_surface_comparison(comparison)

  rows <- lapply(names(comparison$tunings), function(surface_name) {
    tuning <- comparison$tunings[[surface_name]]
    diagnostics <- diagnose_tuning(
      tuning,
      primary_metric = comparison$primary_metric,
      near_best_tolerance = comparison$near_best_tolerance
    )
    overview <- diagnostics$overview
    near_best_fraction <- overview$n_near_best / overview$n_evaluated
    support <- popmaps_tuning_support_label(
      n_near_best = overview$n_near_best,
      near_best_fraction = near_best_fraction
    )

    data.frame(
      surface_name = surface_name,
      validation = overview$validation,
      primary_metric = overview$primary_metric,
      metric_goal = overview$metric_goal,
      n_combinations = overview$n_combinations,
      n_evaluated = overview$n_evaluated,
      n_complete = overview$n_complete,
      best_score = overview$best_score,
      median_score = overview$median_score,
      worst_score = overview$worst_score,
      best_vs_median_percent = overview$best_vs_median_percent,
      best_vs_worst_percent = overview$best_vs_worst_percent,
      near_best_tolerance = overview$near_best_tolerance,
      n_near_best = overview$n_near_best,
      near_best_fraction = near_best_fraction,
      support = support,
      tuning_signal = popmaps_tuning_signal_label(
        best_vs_median_percent = overview$best_vs_median_percent,
        support = support
      ),
      stringsAsFactors = FALSE
    )
  })

  diagnostics <- do.call(rbind, rows)
  rownames(diagnostics) <- NULL
  diagnostics
}

popmaps_surface_near_best_parameter_ranges <- function(comparison) {
  popmaps_check_surface_comparison(comparison)

  rows <- lapply(names(comparison$tunings), function(surface_name) {
    diagnostics <- diagnose_tuning(
      comparison$tunings[[surface_name]],
      primary_metric = comparison$primary_metric,
      near_best_tolerance = comparison$near_best_tolerance
    )
    cbind(
      data.frame(surface_name = surface_name, stringsAsFactors = FALSE),
      diagnostics$parameter_ranges
    )
  })

  ranges <- do.call(rbind, rows)
  rownames(ranges) <- NULL
  ranges
}

popmaps_tuning_support_label <- function(n_near_best, near_best_fraction) {
  if (n_near_best <= 1 || near_best_fraction <= 0.10) {
    return("sharp")
  }
  if (near_best_fraction <= 0.25) {
    return("moderate")
  }

  "broad"
}

popmaps_tuning_signal_label <- function(best_vs_median_percent, support) {
  if (is.finite(best_vs_median_percent) && best_vs_median_percent >= 25 && identical(support, "sharp")) {
    return("strong")
  }
  if ((is.finite(best_vs_median_percent) && best_vs_median_percent >= 10) || !identical(support, "broad")) {
    return("moderate")
  }

  "weak"
}

popmaps_plot_bar <- function(values, colors, ylab, main, ...) {
  if (!any(is.finite(values))) {
    stop("No finite values are available to plot.", call. = FALSE)
  }

  args <- list(...)
  defaults <- list(
    height = values,
    col = colors,
    las = 1,
    ylab = ylab,
    main = main
  )
  defaults[names(args)] <- args
  do.call(graphics::barplot, defaults)
}

popmaps_plot_surface_scores <- function(summary, primary_metric, colors, main, ...) {
  values <- summary$score
  names(values) <- popmaps_surface_plot_labels(summary$surface_name)
  plot_main <- if (is.null(main)) "Surface validation score" else main
  mids <- popmaps_plot_bar(
    values = values,
    colors = colors[summary$surface_name],
    ylab = popmaps_surface_metric_label(primary_metric),
    main = plot_main,
    ...
  )
  graphics::abline(h = values[[1]], lty = 2, col = "gray35")
  graphics::text(
    x = mids,
    y = values,
    labels = paste0("#", summary$rank),
    pos = 3,
    cex = 0.8,
    xpd = NA
  )
}

popmaps_plot_surface_percent_from_best <- function(summary, colors, main, ...) {
  values <- summary$percent_from_best
  names(values) <- popmaps_surface_plot_labels(summary$surface_name)
  plot_main <- if (is.null(main)) "Percent loss from best surface" else main
  mids <- popmaps_plot_bar(
    values = values,
    colors = colors[summary$surface_name],
    ylab = "Percent from best",
    main = plot_main,
    ...
  )
  graphics::abline(h = 5, lty = 2, col = "gray35")
  graphics::text(
    x = mids,
    y = values,
    labels = paste0(signif(values, 3), "%"),
    pos = 3,
    cex = 0.8,
    xpd = NA
  )
}

popmaps_plot_surface_best_parameters <- function(summary, colors, main, ...) {
  parameters <- intersect(
    c("num_sites", "num_tested", "empirical_pt_dist", "popmod"),
    names(summary)
  )
  if (length(parameters) < 1) {
    stop("No best-parameter columns are available to plot.", call. = FALSE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  n_col <- min(2, length(parameters))
  n_row <- ceiling(length(parameters) / n_col)
  graphics::par(mfrow = c(n_row, n_col), mar = c(7, 4, 3, 1), oma = c(0, 0, 3, 0))

  args <- list(...)
  for (parameter in parameters) {
    values <- summary[[parameter]]
    names(values) <- popmaps_surface_plot_labels(summary$surface_name)
    plot_args <- list(
      height = values,
      col = colors[summary$surface_name],
      las = 1,
      ylab = parameter,
      main = parameter
    )
    plot_args[names(args)] <- args
    do.call(graphics::barplot, plot_args)
  }
  if (!is.null(main)) {
    graphics::mtext(main, side = 3, outer = TRUE, line = 1)
  }
}

popmaps_plot_surface_score_distribution <- function(comparison, colors, main, ...) {
  rows <- popmaps_surface_score_distribution(comparison)
  if (!any(is.finite(rows$score))) {
    stop("No finite score values are available to plot.", call. = FALSE)
  }

  surface_levels <- comparison$summary$surface_name
  rows$surface_name <- factor(rows$surface_name, levels = surface_levels)
  labels <- popmaps_surface_plot_labels(surface_levels)
  plot_main <- if (is.null(main)) "Tuning score distribution" else main
  args <- list(...)
  plot_args <- list(
    formula = score ~ surface_name,
    data = rows,
    col = colors[surface_levels],
    las = 2,
    ylab = popmaps_surface_metric_label(comparison$primary_metric),
    xlab = "",
    main = plot_main,
    names = labels
  )
  plot_args[names(args)] <- args
  do.call(graphics::boxplot, plot_args)

  best_scores <- comparison$summary$score[match(surface_levels, comparison$summary$surface_name)]
  graphics::points(seq_along(surface_levels), best_scores, pch = 19, col = colors[surface_levels])
}

popmaps_surface_score_distribution <- function(comparison) {
  rows <- lapply(names(comparison$tunings), function(surface_name) {
    results <- comparison$tunings[[surface_name]]$results
    data.frame(
      surface_name = surface_name,
      score = results[[comparison$primary_metric]],
      failed_folds = results$failed_folds,
      n_scored = results$n_scored,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

popmaps_plot_surface_near_best_parameters <- function(comparison, colors, main, ...) {
  ranges <- popmaps_surface_near_best_parameter_ranges(comparison)
  parameters <- intersect(
    c("num_sites", "num_tested", "empirical_pt_dist", "popmod"),
    unique(ranges$parameter)
  )
  if (length(parameters) < 1) {
    stop("No near-best parameter ranges are available to plot.", call. = FALSE)
  }

  surface_names <- comparison$summary$surface_name
  surface_labels <- popmaps_surface_plot_labels(surface_names)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  n_col <- min(2, length(parameters))
  n_row <- ceiling(length(parameters) / n_col)
  graphics::par(mfrow = c(n_row, n_col), mar = c(7, 4, 3, 1), oma = c(0, 0, 3, 0))

  for (parameter in parameters) {
    parameter_rows <- ranges[ranges$parameter == parameter, , drop = FALSE]
    parameter_rows <- parameter_rows[match(surface_names, parameter_rows$surface_name), , drop = FALSE]
    y_range <- range(
      c(parameter_rows$full_min, parameter_rows$full_max, parameter_rows$near_best_min,
        parameter_rows$near_best_max, parameter_rows$best_value),
      finite = TRUE
    )
    if (!all(is.finite(y_range))) {
      next
    }
    if (diff(y_range) == 0) {
      y_range <- y_range + c(-0.5, 0.5)
    }

    graphics::plot(
      seq_along(surface_names),
      parameter_rows$best_value,
      ylim = y_range,
      xaxt = "n",
      xlab = "",
      ylab = parameter,
      main = parameter,
      pch = 19,
      col = colors[surface_names],
      ...
    )
    graphics::axis(1, at = seq_along(surface_names), labels = surface_labels, las = 1, cex.axis = 0.75)
    graphics::segments(
      seq_along(surface_names),
      parameter_rows$full_min,
      seq_along(surface_names),
      parameter_rows$full_max,
      col = "gray75",
      lwd = 4
    )
    graphics::segments(
      seq_along(surface_names),
      parameter_rows$near_best_min,
      seq_along(surface_names),
      parameter_rows$near_best_max,
      col = colors[surface_names],
      lwd = 3
    )
    graphics::points(
      seq_along(surface_names),
      parameter_rows$best_value,
      pch = 19,
      col = colors[surface_names],
      cex = 1.2
    )
  }

  plot_main <- if (is.null(main)) "Near-best parameter ranges" else main
  graphics::mtext(plot_main, side = 3, outer = TRUE, line = 1)
}

popmaps_check_report_paths <- function(paths, overwrite) {
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0 && !isTRUE(overwrite)) {
    stop(
      "Output files already exist and `overwrite = FALSE`: ",
      paste(existing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

popmaps_write_report_table <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

popmaps_flatten_surface_comparison_grids <- function(grids) {
  rows <- lapply(names(grids), function(surface_name) {
    grid <- grids[[surface_name]]
    data.frame(
      surface_name = surface_name,
      surface = grid$surface,
      surface_values = grid$surface_values,
      distance_units = grid$distance_units,
      reference_distance = grid$reference_distance,
      distance_reference = grid$distance_reference,
      empirical_pt_dist = paste(signif(grid$empirical_pt_dist, 5), collapse = "; "),
      num_sites = paste(grid$num_sites, collapse = "; "),
      num_tested = paste(grid$num_tested, collapse = "; "),
      popmod = paste(signif(grid$popmod, 5), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  grids <- do.call(rbind, rows)
  rownames(grids) <- NULL
  grids
}

popmaps_write_surface_comparison_png <- function(path, width, height, res, expr) {
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(path)
}

popmaps_safe_surface_filename <- function(surface_name) {
  safe <- gsub("[^A-Za-z0-9._-]+", "-", surface_name)
  safe <- gsub("^-+|-+$", "", safe)
  if (!nzchar(safe)) {
    safe <- "surface"
  }

  safe
}

popmaps_write_surface_tuning_results <- function(comparison, dir, prefix, overwrite) {
  paths <- vapply(names(comparison$tunings), function(surface_name) {
    file.path(
      dir,
      paste0(prefix, "-", popmaps_safe_surface_filename(surface_name), "-tuning-results.csv")
    )
  }, character(1))
  popmaps_check_report_paths(paths, overwrite = overwrite)

  for (surface_name in names(comparison$tunings)) {
    rows <- cbind(
      data.frame(surface_name = surface_name, stringsAsFactors = FALSE),
      comparison$tunings[[surface_name]]$results
    )
    popmaps_write_report_table(rows, paths[[surface_name]])
  }

  paths
}

popmaps_markdown_table <- function(x, columns) {
  columns <- intersect(columns, names(x))
  if (nrow(x) < 1 || length(columns) < 1) {
    return("No rows available.")
  }

  x <- x[, columns, drop = FALSE]
  x[] <- lapply(x, function(value) {
    if (is.numeric(value)) {
      signif(value, 5)
    } else {
      value
    }
  })
  x[is.na(x)] <- ""

  header <- paste("|", paste(names(x), collapse = " | "), "|")
  divider <- paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|")
  rows <- apply(x, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  c(header, divider, rows)
}

popmaps_surface_comparison_markdown <- function(comparison,
                                                table_paths,
                                                figure_paths,
                                                tuning_paths) {
  metric_sentence <- if (popmaps_metric_is_maximized(comparison$primary_metric)) {
    "Higher values are better for the primary metric."
  } else {
    "Lower values are better for the primary metric."
  }

  report <- c(
    "# POPMAPS Surface Comparison Report",
    "",
    paste0("Decision: `", comparison$support$decision, "`"),
    paste0("Best surface: `", comparison$support$best_surface, "`"),
    paste0("Primary metric: `", comparison$primary_metric, "`"),
    paste0("Validation: `", comparison$validation, "`"),
    "",
    metric_sentence,
    "`percent_from_best` is the relative change from the best-ranked surface. ",
    "Surfaces inside the near-best tolerance should be treated as similarly supported by the current validation design.",
    "",
    "## Tables",
    "",
    paste0("- Summary: `", basename(table_paths[["summary"]]), "`"),
    paste0("- Support: `", basename(table_paths[["support"]]), "`"),
    paste0("- Tuning grids: `", basename(table_paths[["grids"]]), "`"),
    paste0("- Tuning diagnostics: `", basename(table_paths[["tuning_diagnostics"]]), "`"),
    paste0("- Near-best parameter ranges: `", basename(table_paths[["near_best_parameters"]]), "`"),
    "",
    "## Surface Support",
    "",
    popmaps_markdown_table(
      comparison$support,
      c("decision", "best_surface", "primary_metric", "best_score",
        "near_best_tolerance", "n_near_best", "near_best_surfaces")
    ),
    "",
    "## Surface Ranking",
    "",
    popmaps_markdown_table(
      comparison$summary,
      c("rank", "surface_name", "surface", "surface_values", "score",
        "percent_from_best", "distance_units", "num_sites", "num_tested",
        "popmod", "empirical_pt_dist", "failed_folds")
    ),
    "",
    "## Tuning Diagnostics",
    "",
    popmaps_markdown_table(
      popmaps_surface_tuning_diagnostics(comparison),
      c("surface_name", "support", "tuning_signal", "n_near_best",
        "near_best_fraction", "best_score", "median_score",
        "best_vs_median_percent")
    ),
    "",
    "## Figures",
    "",
    paste0("- `figures/", basename(figure_paths[["score"]]), "`: primary validation score by candidate surface."),
    paste0("- `figures/", basename(figure_paths[["percent_from_best"]]), "`: relative loss from the best surface."),
    paste0("- `figures/", basename(figure_paths[["best_parameters"]]), "`: best parameter values selected for each surface."),
    paste0("- `figures/", basename(figure_paths[["score_distribution"]]), "`: tuning score distribution across parameter combinations."),
    paste0("- `figures/", basename(figure_paths[["near_best_parameters"]]), "`: full and near-best parameter ranges by surface.")
  )

  if (length(tuning_paths) > 0) {
    report <- c(
      report,
      "",
      "## Full Tuning Results",
      "",
      paste0("- `", basename(tuning_paths), "`")
    )
  }

  report
}
