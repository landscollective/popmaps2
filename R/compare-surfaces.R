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
      surface_grid = surface_grid
    )
    tunings[[surface_name]] <- tune_popmaps(
      input_raster = surface_object,
      input_locs = input_locs,
      surface = surface_object$surface,
      empirical_pt_dist = grids[[surface_name]]$empirical_pt_dist,
      num_sites = grids[[surface_name]]$num_sites,
      num_tested = grids[[surface_name]]$num_tested,
      popmod = grids[[surface_name]]$popmod,
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
    grids = grids,
    primary_metric = primary_metric,
    validation = validation,
    surface_grid = surface_grid,
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

popmaps_resolve_comparison_grid <- function(input_locs,
                                            surface_object,
                                            num_sites,
                                            num_tested,
                                            popmod,
                                            empirical_pt_dist,
                                            surface_grid) {
  suggested <- if (surface_grid == "surface_specific") {
    suggest_surface_tuning_grid(
      input_raster = surface_object,
      input_locs = input_locs,
      surface = surface_object$surface,
      surface_values = surface_object$surface_values,
      num_sites = num_sites,
      num_tested = num_tested
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
#'   surface.
#' @param col Optional vector of plotting colors. Named vectors are matched to
#'   surface names; unnamed vectors are recycled in comparison order.
#' @param main Optional plot title.
#' @param ... Additional arguments passed to base plotting functions.
#'
#' @return Invisibly returns `comparison`.
#'
#' @export
plot_surface_comparison <- function(comparison,
                                    type = c("score", "percent_from_best", "best_parameters"),
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
  } else {
    popmaps_plot_surface_best_parameters(summary, colors, main, ...)
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
    grids = file.path(dir, paste0(prefix, "-grids.csv"))
  )
  figure_paths <- c(
    score = file.path(figure_dir, paste0(prefix, "-scores.png")),
    percent_from_best = file.path(figure_dir, paste0(prefix, "-percent-from-best.png")),
    best_parameters = file.path(figure_dir, paste0(prefix, "-best-parameters.png"))
  )
  report_path <- file.path(dir, paste0(prefix, "-report.md"))

  popmaps_check_report_paths(
    c(table_paths, figure_paths, report = report_path),
    overwrite = overwrite
  )

  grids <- popmaps_flatten_surface_comparison_grids(comparison$grids)
  popmaps_write_report_table(comparison$summary, table_paths[["summary"]])
  popmaps_write_report_table(comparison$support, table_paths[["support"]])
  popmaps_write_report_table(grids, table_paths[["grids"]])

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

popmaps_plot_bar <- function(values, colors, ylab, main, ...) {
  if (!any(is.finite(values))) {
    stop("No finite values are available to plot.", call. = FALSE)
  }

  args <- list(...)
  defaults <- list(
    height = values,
    col = colors,
    las = 2,
    ylab = ylab,
    main = main
  )
  defaults[names(args)] <- args
  do.call(graphics::barplot, defaults)
}

popmaps_plot_surface_scores <- function(summary, primary_metric, colors, main, ...) {
  values <- summary$score
  names(values) <- summary$surface_name
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
  names(values) <- summary$surface_name
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
    names(values) <- summary$surface_name
    plot_args <- list(
      height = values,
      col = colors[summary$surface_name],
      las = 2,
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
    "## Figures",
    "",
    paste0("- `figures/", basename(figure_paths[["score"]]), "`: primary validation score by candidate surface."),
    paste0("- `figures/", basename(figure_paths[["percent_from_best"]]), "`: relative loss from the best surface."),
    paste0("- `figures/", basename(figure_paths[["best_parameters"]]), "`: best parameter values selected for each surface.")
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
