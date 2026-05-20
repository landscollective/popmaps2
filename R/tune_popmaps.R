#' Tune POPMAPS parameters with cross-validation
#'
#' @description `tune_popmaps()` compares POPMAPS parameter combinations by
#'   withholding empirical sampling sites, predicting ancestry coefficients from
#'   the remaining sites, and summarizing prediction error. It supports
#'   geographic-distance tuning and modern least-cost tuning for prepared
#'   suitability, conductance, or resistance surfaces.
#'
#' @param input_raster A `terra::SpatRaster`, `raster::RasterLayer`, or path to
#'   a raster file defining the interpolation surface. It may also be a
#'   `popmaps_surface` object returned by [prepare_popmaps_surface()].
#' @param input_locs A data frame or matrix with sampling location name,
#'   longitude, latitude, and one or more ancestry coefficient columns.
#' @param surface Character. `"G"` tunes geographic-distance interpolation.
#'   `"C"` tunes suitability- or conductance-weighted least-cost interpolation
#'   using the modern internal distance helper.
#' @param surface_values Character. Meaning of raster values when
#'   `surface = "C"`. `"suitability"` and `"conductance"` use values directly;
#'   `"resistance"` converts values to conductance before distances are
#'   calculated. Ignored when `input_raster` is already a `popmaps_surface`
#'   object.
#' @param empirical_pt_dist Numeric vector. Minimum distances required between
#'   empirical sites selected for a prediction. Values are kilometers for
#'   `surface = "G"` and least-cost distance units for `surface = "C"`.
#' @param num_sites Integer vector. Candidate pool sizes to evaluate.
#' @param num_tested Integer vector. Numbers of empirical sites used to estimate
#'   ancestry coefficients.
#' @param popmod Numeric vector. Distance-decay parameter values to evaluate.
#' @param threshold Numeric scalar. Raster values below this threshold are not
#'   scored.
#' @param validation Cross-validation design. `"loo"` withholds one site at a
#'   time. `"spatial_block"` withholds spatially grouped sites, which is a
#'   stricter test of prediction into undersampled regions.
#' @param n_blocks Target number of spatial blocks when
#'   `validation = "spatial_block"` and `block_assignments = NULL`.
#' @param block_assignments Optional vector assigning each empirical site to a
#'   spatial block. If supplied, it must have one value per row in `input_locs`.
#' @param spatial_block_repeats Number of spatial-block layouts to evaluate when
#'   `validation = "spatial_block"` and `block_assignments = NULL`. Values
#'   greater than one repeat the spatial-block validation with rotated spatial
#'   partitions and report repeat-level uncertainty.
#' @param spatial_block_seed Optional random seed for repeated spatial-block
#'   layouts. The first repeat uses the deterministic default partition; later
#'   repeats use random spatial rotations.
#' @param primary_metric Metric used to select the best parameter combination.
#' @param dist_prob_func Function defining the relationship between distance and
#'   empirical-site contribution.
#' @param rescale_conductance Logical. If `TRUE`, scale conductance values by
#'   the largest non-missing conductance value before least-cost distances are
#'   calculated for `surface = "C"`.
#' @param resistance_epsilon Positive numeric scalar added to resistance values
#'   before inversion when `surface_values = "resistance"`.
#' @param quiet Logical. If `FALSE`, print a short completion message.
#'
#' @return A `popmaps_tuning` object with:
#' \describe{
#'   \item{results}{One row per parameter combination with average validation metrics.}
#'   \item{folds}{One row per withheld site and parameter combination.}
#'   \item{best}{The best row from `results` according to `primary_metric`.}
#' }
#'
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of
#'   genetic diversity using ancestry probability surfaces. Methods in Ecology
#'   and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#'
#' @examples
#' ex_raster <- raster::aggregate(hija_raster, fact = 64)
#' tuning <- tune_popmaps(
#'   input_raster = ex_raster,
#'   input_locs = hija_struc,
#'   num_sites = c(8, 10),
#'   num_tested = c(3, 4),
#'   popmod = c(-0.01, -0.05),
#'   quiet = TRUE
#' )
#' tuning$best
#'
#' @export
tune_popmaps <- function(input_raster = "",
                         input_locs = "",
                         surface = "G",
                         surface_values = c("suitability", "conductance", "resistance"),
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
                         rescale_conductance = FALSE,
                         resistance_epsilon = sqrt(.Machine$double.eps),
                         quiet = TRUE) {
  surface <- match.arg(surface, c("G", "C"))
  surface_values <- match.arg(surface_values)
  validation <- match.arg(validation)
  primary_metric <- match.arg(primary_metric)

  if (!is.function(dist_prob_func)) {
    stop("`dist_prob_func` must be a function.", call. = FALSE)
  }
  if (!is.logical(rescale_conductance) || length(rescale_conductance) != 1 || is.na(rescale_conductance)) {
    stop("`rescale_conductance` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  num_sites <- popmaps_check_tuning_values(
    num_sites,
    "`num_sites`",
    whole_number = TRUE,
    positive = TRUE
  )
  num_tested <- popmaps_check_tuning_values(
    num_tested,
    "`num_tested`",
    whole_number = TRUE,
    positive = TRUE
  )
  popmod <- popmaps_check_tuning_values(popmod, "`popmod`")
  popmaps_check_finite_scalar(threshold, "`threshold`")
  empirical_pt_dist <- popmaps_check_tuning_values(
    empirical_pt_dist,
    "`empirical_pt_dist`",
    nonnegative = TRUE
  )
  popmaps_check_whole_count(spatial_block_repeats, "`spatial_block_repeats`", positive = TRUE)
  if (!is.null(spatial_block_seed)) {
    popmaps_check_whole_count(spatial_block_seed, "`spatial_block_seed`", allow_zero = TRUE)
  }

  parameter_grid <- popmaps_make_tuning_grid(
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    empirical_pt_dist = empirical_pt_dist
  )

  popmaps_evaluate_tuning_grid(
    input_raster = input_raster,
    input_locs = input_locs,
    surface = surface,
    parameter_grid = parameter_grid,
    threshold = threshold,
    validation = validation,
    n_blocks = n_blocks,
    block_assignments = block_assignments,
    spatial_block_repeats = spatial_block_repeats,
    spatial_block_seed = spatial_block_seed,
    primary_metric = primary_metric,
    dist_prob_func = dist_prob_func,
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon,
    quiet = quiet,
    call = match.call()
  )
}

#' Suggest a data-adaptive POPMAPS tuning grid
#'
#' @description `suggest_tuning_grid()` uses empirical sampling-site distances
#'   to suggest candidate values for `tune_popmaps()` and
#'   `adaptive_tune_popmaps()`. The suggestions are intended as a transparent
#'   starting point, not as universal defaults.
#'
#' @param input_locs A data frame or matrix with sampling location name,
#'   longitude, latitude, and one or more ancestry coefficient columns.
#' @param num_sites Optional integer vector of candidate site-pool sizes.
#' @param num_tested Optional integer vector of candidate numbers of sites used
#'   in each prediction.
#' @param empirical_pt_dist_probs Distance quantiles used to suggest
#'   `empirical_pt_dist`, after always including zero.
#' @param distance_weights Desired weights retained at the reference distance;
#'   converted to `popmod` values with `log(weight) / reference_distance`.
#' @param distance_reference Which site-distance summary to use as the
#'   reference distance for `popmod`.
#' @param max_num_tested Maximum automatically suggested `num_tested` value.
#'
#' @return A `popmaps_tuning_grid` list with suggested parameter vectors and
#'   empirical distance summaries.
#'
#' @examples
#' grid <- suggest_tuning_grid(hija_struc)
#' grid
#'
#' @export
suggest_tuning_grid <- function(input_locs,
                                num_sites = NULL,
                                num_tested = NULL,
                                empirical_pt_dist_probs = c(0.05, 0.10, 0.25),
                                distance_weights = c(0.95, 0.75, 0.50, 0.25, 0.10, 0.05),
                                distance_reference = c("median", "mean"),
                                max_num_tested = 8) {
  distance_reference <- match.arg(distance_reference)
  locations <- popmaps_prepare_locations(input_locs)
  coords <- as.matrix(locations[, 2:3, drop = FALSE])
  distance_matrix <- popmaps_empirical_site_distances(coords)

  popmaps_suggest_tuning_grid_from_distances(
    locations = locations,
    distance_matrix = distance_matrix,
    num_sites = num_sites,
    num_tested = num_tested,
    empirical_pt_dist_probs = empirical_pt_dist_probs,
    distance_weights = distance_weights,
    distance_reference = distance_reference,
    max_num_tested = max_num_tested,
    surface = "G",
    surface_values = NA_character_,
    distance_units = "km"
  )
}

#' Suggest a POPMAPS tuning grid from a candidate surface
#'
#' @description
#' `suggest_surface_tuning_grid()` is the surface-aware companion to
#' [suggest_tuning_grid()]. It derives `empirical_pt_dist` and `popmod` from
#' empirical sampling-site distances measured over the selected surface. For
#' `surface = "G"` those distances are geographic kilometers. For
#' `surface = "C"` they are least-cost distances through a suitability,
#' conductance, or resistance surface.
#'
#' @inheritParams tune_popmaps
#' @inheritParams suggest_tuning_grid
#'
#' @return A `popmaps_tuning_grid` list.
#'
#' @examples
#' ex_raster <- terra::rast(raster::aggregate(hija_raster, fact = 240))
#' grid <- suggest_surface_tuning_grid(
#'   input_raster = ex_raster,
#'   input_locs = hija_struc,
#'   surface = "G"
#' )
#' grid$distance_units
#'
#' @export
suggest_surface_tuning_grid <- function(input_raster,
                                        input_locs,
                                        surface = "G",
                                        surface_values = c("suitability", "conductance", "resistance"),
                                        num_sites = NULL,
                                        num_tested = NULL,
                                        empirical_pt_dist_probs = c(0.05, 0.10, 0.25),
                                        distance_weights = c(0.95, 0.75, 0.50, 0.25, 0.10, 0.05),
                                        distance_reference = c("median", "mean"),
                                        max_num_tested = 8,
                                        rescale_conductance = FALSE,
                                        resistance_epsilon = sqrt(.Machine$double.eps)) {
  surface <- match.arg(surface, c("G", "C"))
  surface_values <- if (surface == "C") match.arg(surface_values) else NA_character_
  distance_reference <- match.arg(distance_reference)
  if (!is.logical(rescale_conductance) || length(rescale_conductance) != 1 || is.na(rescale_conductance)) {
    stop("`rescale_conductance` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  locations <- popmaps_prepare_locations(input_locs)
  coords <- as.matrix(locations[, 2:3, drop = FALSE])

  if (surface == "G") {
    distance_matrix <- popmaps_empirical_site_distances(coords)
    grid_surface_values <- NA_character_
    distance_units <- "km"
  } else {
    surface_object <- if (inherits(input_raster, "popmaps_surface")) {
      if (!identical(input_raster$surface, "C")) {
        stop("A prepared surface object must have `surface = \"C\"` for least-cost grid suggestions.", call. = FALSE)
      }
      input_raster
    } else {
      prepare_popmaps_surface(
        input_raster = input_raster,
        surface = "C",
        surface_values = surface_values,
        rescale_conductance = rescale_conductance,
        resistance_epsilon = resistance_epsilon
      )
    }
    graph <- popmaps_cost_distance_graph(surface_object, directions = 8)
    distance_matrix <- popmaps_cost_distance_matrix(
      surface = surface_object,
      from_coords = coords,
      directions = 8,
      graph = graph
    )
    grid_surface_values <- surface_object$surface_values
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
    surface = surface,
    surface_values = grid_surface_values,
    distance_units = distance_units
  )
}

popmaps_suggest_tuning_grid_from_distances <- function(locations,
                                                       distance_matrix,
                                                       num_sites,
                                                       num_tested,
                                                       empirical_pt_dist_probs,
                                                       distance_weights,
                                                       distance_reference,
                                                       max_num_tested,
                                                       surface,
                                                       surface_values,
                                                       distance_units) {
  if (nrow(locations) < 3) {
    stop("At least three empirical sites are required for leave-one-site-out tuning.", call. = FALSE)
  }
  if (!is.matrix(distance_matrix) || nrow(distance_matrix) != nrow(locations) || ncol(distance_matrix) != nrow(locations)) {
    stop("`distance_matrix` must be a square matrix with one row and column per empirical site.", call. = FALSE)
  }
  distances <- distance_matrix[upper.tri(distance_matrix)]
  distances <- distances[is.finite(distances) & distances > 0]
  if (length(distances) < 1) {
    stop("Empirical sites must include at least two distinct locations with finite positive distances.", call. = FALSE)
  }

  n_training <- nrow(locations) - 1
  if (is.null(num_sites)) {
    num_sites <- popmaps_default_num_sites(n_training)
  } else {
    num_sites <- popmaps_check_tuning_values(
      num_sites,
      "`num_sites`",
      whole_number = TRUE,
      positive = TRUE
    )
    if (max(num_sites) > n_training) {
      stop("`num_sites` cannot exceed the number of leave-one-out training sites.", call. = FALSE)
    }
  }

  if (is.null(num_tested)) {
    popmaps_check_whole_count(max_num_tested, "`max_num_tested`", positive = TRUE)
    max_tested <- min(max_num_tested, min(num_sites))
    if (max_tested < 2) {
      stop("Suggested grids require `num_sites` to allow at least two tested sites.", call. = FALSE)
    }
    num_tested <- seq.int(2, max_tested)
  } else {
    num_tested <- popmaps_check_tuning_values(
      num_tested,
      "`num_tested`",
      whole_number = TRUE,
      positive = TRUE
    )
  }
  if (max(num_tested) > min(num_sites)) {
    stop("Suggested grids require all `num_tested` values to be <= all `num_sites` values.", call. = FALSE)
  }

  if (
    !is.numeric(empirical_pt_dist_probs) ||
      length(empirical_pt_dist_probs) < 1 ||
      any(!is.finite(empirical_pt_dist_probs)) ||
      any(empirical_pt_dist_probs < 0 | empirical_pt_dist_probs > 1)
  ) {
    stop("`empirical_pt_dist_probs` must contain probabilities between 0 and 1.", call. = FALSE)
  }
  empirical_pt_dist <- sort(unique(c(
    0,
    as.numeric(stats::quantile(
      distances,
      probs = empirical_pt_dist_probs,
      names = FALSE,
      type = 7
    ))
  )))

  if (
    !is.numeric(distance_weights) ||
      length(distance_weights) < 1 ||
      any(!is.finite(distance_weights)) ||
      any(distance_weights <= 0 | distance_weights >= 1)
  ) {
    stop("`distance_weights` must contain finite values between 0 and 1.", call. = FALSE)
  }
  reference_distance <- if (distance_reference == "median") {
    stats::median(distances)
  } else {
    mean(distances)
  }
  if (!is.finite(reference_distance) || reference_distance <= 0) {
    stop("The empirical reference distance must be positive.", call. = FALSE)
  }
  popmod <- sort(unique(as.numeric(log(distance_weights) / reference_distance)))

  grid <- list(
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    empirical_pt_dist = empirical_pt_dist,
    distance_summary = stats::quantile(
      distances,
      probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1),
      names = TRUE,
      type = 7
    ),
    distance_reference = distance_reference,
    reference_distance = reference_distance,
    distance_weights = distance_weights,
    surface = surface,
    surface_values = surface_values,
    distance_units = distance_units
  )
  class(grid) <- "popmaps_tuning_grid"
  grid
}

#' Adaptively tune POPMAPS parameters
#'
#' @description `adaptive_tune_popmaps()` samples parameter combinations from a
#'   data-adaptive search space, evaluates them with the same leave-one-site-out
#'   scoring used by `tune_popmaps()`, then optionally refines sampling around
#'   the best-performing region.
#'
#' @param input_raster A `terra::SpatRaster`, `raster::RasterLayer`, or path to
#'   a raster file defining the interpolation surface.
#' @param input_locs A data frame or matrix with sampling location name,
#'   longitude, latitude, and one or more ancestry coefficient columns.
#' @param surface Character. Currently only `"G"` is supported.
#' @param empirical_pt_dist Optional numeric vector or range for the rarefaction
#'   distance search space. If `NULL`, values are suggested from site distances.
#' @param num_sites Optional integer candidate site-pool sizes.
#' @param num_tested Optional integer candidate numbers of sites used in each
#'   prediction.
#' @param popmod Optional numeric vector or range for the distance-decay search
#'   space. If `NULL`, values are suggested from site distances.
#' @param n_initial Number of initial parameter combinations to sample.
#' @param n_refine Number of additional combinations to sample around the best
#'   initial region.
#' @param method Search method. `"latin_hypercube"` stratifies samples across
#'   each parameter range; `"random"` samples independently.
#' @param seed Optional random seed for reproducibility.
#' @inheritParams tune_popmaps
#'
#' @return A `popmaps_adaptive_tuning` object. It has the same `results`,
#'   `folds`, and `best` elements as `tune_popmaps()`, plus a `search` element
#'   describing the sampled search.
#'
#' @examples
#' ex_raster <- raster::aggregate(hija_raster, fact = 240)
#' adaptive <- adaptive_tune_popmaps(
#'   input_raster = ex_raster,
#'   input_locs = hija_struc,
#'   n_initial = 6,
#'   n_refine = 4,
#'   seed = 1,
#'   quiet = TRUE
#' )
#' adaptive$best
#'
#' @export
adaptive_tune_popmaps <- function(input_raster = "",
                                  input_locs = "",
                                  surface = "G",
                                  empirical_pt_dist = NULL,
                                  num_sites = NULL,
                                  num_tested = NULL,
                                  popmod = NULL,
                                  n_initial = 50,
                                  n_refine = 50,
                                  method = c("latin_hypercube", "random"),
                                  seed = NULL,
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
                                  quiet = TRUE) {
  surface <- match.arg(surface, c("G", "C"))
  method <- match.arg(method)
  validation <- match.arg(validation)
  primary_metric <- match.arg(primary_metric)
  if (surface != "G") {
    stop(
      "`adaptive_tune_popmaps()` currently supports only geographic-distance tuning with surface = 'G'.",
      call. = FALSE
    )
  }
  if (!is.function(dist_prob_func)) {
    stop("`dist_prob_func` must be a function.", call. = FALSE)
  }
  popmaps_check_finite_scalar(threshold, "`threshold`")
  popmaps_check_whole_count(n_initial, "`n_initial`", positive = TRUE)
  popmaps_check_whole_count(n_refine, "`n_refine`", allow_zero = TRUE)
  popmaps_check_whole_count(spatial_block_repeats, "`spatial_block_repeats`", positive = TRUE)
  if (!is.null(spatial_block_seed)) {
    popmaps_check_whole_count(spatial_block_seed, "`spatial_block_seed`", allow_zero = TRUE)
  }
  if (!is.null(seed)) {
    popmaps_check_whole_count(seed, "`seed`", allow_zero = TRUE)
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  suggested <- suggest_tuning_grid(
    input_locs = input_locs
  )
  search_space <- list(
    num_sites = popmaps_resolve_search_values(num_sites, suggested$num_sites, "`num_sites`",
      whole_number = TRUE,
      positive = TRUE
    ),
    num_tested = popmaps_resolve_search_values(num_tested, suggested$num_tested, "`num_tested`",
      whole_number = TRUE,
      positive = TRUE
    ),
    popmod = popmaps_resolve_search_values(popmod, suggested$popmod, "`popmod`"),
    empirical_pt_dist = popmaps_resolve_search_values(
      empirical_pt_dist,
      suggested$empirical_pt_dist,
      "`empirical_pt_dist`",
      nonnegative = TRUE
    )
  )

  initial_grid <- popmaps_sample_tuning_space(
    n = n_initial,
    search_space = search_space,
    method = method
  )
  initial_tuning <- popmaps_evaluate_tuning_grid(
    input_raster = input_raster,
    input_locs = input_locs,
    surface = surface,
    parameter_grid = initial_grid,
    threshold = threshold,
    validation = validation,
    n_blocks = n_blocks,
    block_assignments = block_assignments,
    spatial_block_repeats = spatial_block_repeats,
    spatial_block_seed = spatial_block_seed,
    primary_metric = primary_metric,
    dist_prob_func = dist_prob_func,
    quiet = TRUE,
    call = match.call()
  )

  if (n_refine > 0) {
    refine_space <- popmaps_refine_tuning_space(
      results = initial_tuning$results,
      search_space = search_space,
      primary_metric = primary_metric
    )
    refine_grid <- popmaps_sample_tuning_space(
      n = n_refine,
      search_space = refine_space,
      method = method
    )
    parameter_grid <- popmaps_normalize_tuning_grid(rbind(
      initial_grid[, c("num_sites", "num_tested", "popmod", "empirical_pt_dist")],
      refine_grid[, c("num_sites", "num_tested", "popmod", "empirical_pt_dist")]
    ))
  } else {
    refine_grid <- NULL
    parameter_grid <- initial_grid
  }

  search <- list(
    method = method,
    n_initial = n_initial,
    n_refine = n_refine,
    initial_candidates = nrow(initial_grid),
    refine_candidates = if (is.null(refine_grid)) 0L else nrow(refine_grid),
    final_candidates = nrow(parameter_grid),
    suggested_grid = suggested
  )

  tuning <- popmaps_evaluate_tuning_grid(
    input_raster = input_raster,
    input_locs = input_locs,
    surface = surface,
    parameter_grid = parameter_grid,
    threshold = threshold,
    validation = validation,
    n_blocks = n_blocks,
    block_assignments = block_assignments,
    spatial_block_repeats = spatial_block_repeats,
    spatial_block_seed = spatial_block_seed,
    primary_metric = primary_metric,
    dist_prob_func = dist_prob_func,
    quiet = quiet,
    call = match.call(),
    class = c("popmaps_adaptive_tuning", "popmaps_tuning"),
    search = search
  )

  tuning
}

#' Diagnose a POPMAPS tuning result
#'
#' @description `diagnose_tuning()` summarizes whether parameter tuning found a
#'   clearly supported parameter combination or a broad set of similarly
#'   performing alternatives. It is intended to help users interpret tuning
#'   output biologically instead of choosing a row from `tuning$results` by eye.
#'
#' @param tuning A `popmaps_tuning` or `popmaps_adaptive_tuning` object returned
#'   by `tune_popmaps()` or `adaptive_tune_popmaps()`.
#' @param primary_metric Metric used to rank parameter combinations. Defaults to
#'   the metric stored in `tuning`.
#' @param near_best_tolerance Non-negative relative tolerance used to define
#'   near-best parameter combinations. The default, `0.05`, keeps combinations
#'   within 5% of the best score for the primary metric.
#' @param complete_only Logical. If `TRUE`, diagnose only parameter combinations
#'   with no failed validation folds.
#'
#' @return A `popmaps_tuning_diagnostics` list with:
#' \describe{
#'   \item{overview}{One-row summary of tuning strength and near-best support.}
#'   \item{near_best}{Parameter combinations within `near_best_tolerance` of the best score.}
#'   \item{parameter_ranges}{Near-best and full-grid support for each tuning parameter.}
#'   \item{parameter_effects}{Average score by parameter value.}
#' }
#'
#' @examples
#' ex_raster <- raster::aggregate(hija_raster, fact = 240)
#' tuning <- tune_popmaps(
#'   input_raster = ex_raster,
#'   input_locs = hija_struc,
#'   empirical_pt_dist = c(0, 5),
#'   num_sites = c(5, 6),
#'   num_tested = c(2, 3),
#'   popmod = c(-0.01, -0.05),
#'   quiet = TRUE
#' )
#' diagnose_tuning(tuning)
#'
#' @export
diagnose_tuning <- function(tuning,
                            primary_metric = tuning$primary_metric,
                            near_best_tolerance = 0.05,
                            complete_only = TRUE) {
  if (!is.list(tuning) || is.null(tuning$results)) {
    stop("`tuning` must be a tuning object returned by `tune_popmaps()`.", call. = FALSE)
  }
  if (!is.character(primary_metric) || length(primary_metric) != 1) {
    stop("`primary_metric` must be one metric name.", call. = FALSE)
  }
  popmaps_check_finite_scalar(near_best_tolerance, "`near_best_tolerance`")
  if (near_best_tolerance < 0) {
    stop("`near_best_tolerance` must be non-negative.", call. = FALSE)
  }
  if (!is.logical(complete_only) || length(complete_only) != 1 || is.na(complete_only)) {
    stop("`complete_only` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  results <- tuning$results
  if (!primary_metric %in% names(results)) {
    stop("`primary_metric` must be a column in `tuning$results`.", call. = FALSE)
  }

  if (isTRUE(complete_only) && all(c("failed_folds", "n_scored") %in% names(results))) {
    results <- results[results$failed_folds == 0 & results$n_scored > 0, , drop = FALSE]
  }
  if (nrow(results) < 1) {
    stop("No complete tuning results are available to diagnose.", call. = FALSE)
  }

  metric_values <- results[[primary_metric]]
  scored <- is.finite(metric_values)
  if (!any(scored)) {
    stop("No finite values are available for `primary_metric`.", call. = FALSE)
  }

  maximize <- popmaps_metric_is_maximized(primary_metric)
  best <- popmaps_tuning_best_row(results, primary_metric)
  best_score <- best[[primary_metric]][1]
  metric_scale <- max(abs(best_score), .Machine$double.eps)
  near_best <- if (maximize) {
    scored & metric_values >= best_score - near_best_tolerance * metric_scale
  } else {
    scored & metric_values <= best_score + near_best_tolerance * metric_scale
  }
  near_best_results <- results[near_best, , drop = FALSE]
  near_best_results <- near_best_results[
    popmaps_order_tuning_results(near_best_results, primary_metric),
    ,
    drop = FALSE
  ]

  median_score <- stats::median(metric_values[scored])
  worst_score <- if (maximize) min(metric_values[scored]) else max(metric_values[scored])
  best_vs_median_delta <- if (maximize) best_score - median_score else median_score - best_score
  best_vs_worst_delta <- if (maximize) best_score - worst_score else worst_score - best_score

  overview <- data.frame(
    validation = if (!is.null(tuning$validation)) tuning$validation else paste(unique(results$validation), collapse = ", "),
    primary_metric = primary_metric,
    metric_goal = if (maximize) "maximize" else "minimize",
    n_combinations = nrow(tuning$results),
    n_evaluated = nrow(results),
    n_complete = if ("failed_folds" %in% names(results)) sum(results$failed_folds == 0) else NA_integer_,
    best_score = best_score,
    median_score = median_score,
    worst_score = worst_score,
    best_vs_median_delta = best_vs_median_delta,
    best_vs_median_percent = popmaps_percent_change(best_vs_median_delta, median_score),
    best_vs_worst_delta = best_vs_worst_delta,
    best_vs_worst_percent = popmaps_percent_change(best_vs_worst_delta, worst_score),
    near_best_tolerance = near_best_tolerance,
    n_near_best = nrow(near_best_results),
    stringsAsFactors = FALSE
  )

  tuning_parameters <- intersect(
    c("num_sites", "num_tested", "popmod", "half_distance",
      "ten_pct_distance", "empirical_pt_dist"),
    names(results)
  )

  parameter_ranges <- do.call(rbind, lapply(tuning_parameters, function(parameter) {
    full_values <- results[[parameter]]
    near_values <- near_best_results[[parameter]]

    data.frame(
      parameter = parameter,
      best_value = best[[parameter]][1],
      near_best_min = min(near_values, na.rm = TRUE),
      near_best_max = max(near_values, na.rm = TRUE),
      near_best_unique = paste(signif(sort(unique(near_values)), 5), collapse = ", "),
      full_min = min(full_values, na.rm = TRUE),
      full_max = max(full_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(parameter_ranges) <- NULL

  effect_parameters <- intersect(
    c("num_sites", "num_tested", "popmod", "empirical_pt_dist"),
    names(results)
  )
  parameter_effects <- do.call(rbind, lapply(effect_parameters, function(parameter) {
    values <- sort(unique(results[[parameter]]))
    rows <- lapply(values, function(value) {
      idx <- results[[parameter]] == value
      scores <- results[[primary_metric]][idx]
      finite_scores <- scores[is.finite(scores)]
      mean_score <- mean(scores, na.rm = TRUE)
      if (is.nan(mean_score)) {
        mean_score <- NA_real_
      }
      score_sd <- if (length(finite_scores) > 1) stats::sd(finite_scores) else NA_real_

      data.frame(
        parameter = parameter,
        value = value,
        n_combinations = sum(idx),
        n_scored = sum(is.finite(scores)),
        mean_score = mean_score,
        score_sd = score_sd,
        loss_from_best = if (maximize) best_score - mean_score else mean_score - best_score,
        is_best_value = value == best[[parameter]][1],
        stringsAsFactors = FALSE
      )
    })
    parameter_rows <- do.call(rbind, rows)
    parameter_rows <- parameter_rows[
      order(if (maximize) -parameter_rows$mean_score else parameter_rows$mean_score),
      ,
      drop = FALSE
    ]
    parameter_rows$rank <- seq_len(nrow(parameter_rows))
    parameter_rows
  }))
  rownames(parameter_effects) <- NULL

  diagnostics <- list(
    overview = overview,
    near_best = near_best_results,
    parameter_ranges = parameter_ranges,
    parameter_effects = parameter_effects
  )
  class(diagnostics) <- "popmaps_tuning_diagnostics"
  diagnostics
}

#' @export
print.popmaps_tuning_diagnostics <- function(x, ...) {
  cat("POPMAPS tuning diagnostics\n")
  cat("Primary metric: ", x$overview$primary_metric, " (", x$overview$metric_goal, ")\n", sep = "")
  cat("Validation: ", x$overview$validation, "\n", sep = "")
  cat("Near-best combinations: ", x$overview$n_near_best, "\n\n", sep = "")
  cat("Overview:\n")
  print(x$overview, row.names = FALSE)
  cat("\nNear-best parameter support:\n")
  print(x$parameter_ranges, row.names = FALSE)
  invisible(x)
}

popmaps_evaluate_tuning_grid <- function(input_raster,
                                         input_locs,
                                         surface,
                                         parameter_grid,
                                         threshold,
                                         validation,
                                         n_blocks,
                                         block_assignments,
                                         spatial_block_repeats,
                                         spatial_block_seed,
                                         primary_metric,
                                         dist_prob_func,
                                         surface_values = "suitability",
                                         rescale_conductance = FALSE,
                                         resistance_epsilon = sqrt(.Machine$double.eps),
                                         quiet,
                                         call,
                                         class = "popmaps_tuning",
                                         search = NULL) {
  if (inherits(input_raster, "popmaps_surface")) {
    surface_object <- input_raster
    if (!identical(surface_object$surface, surface)) {
      stop("`input_raster` surface metadata does not match `surface`.", call. = FALSE)
    }
    prepare_raster <- surface_object$rast
  } else {
    surface_object <- NULL
    prepare_raster <- input_raster
  }

  prepared <- popmaps_prepare_inputs(
    input_raster = prepare_raster,
    input_locs = input_locs,
    surface = surface,
    num_sites = max(parameter_grid$num_sites),
    num_tested = max(parameter_grid$num_tested),
    threshold = threshold,
    empirical_pt_dist = max(parameter_grid$empirical_pt_dist),
    jackknife = TRUE,
    require_legacy_c = FALSE
  )

  locations <- prepared$locations
  coords <- as.matrix(locations[, 2:3, drop = FALSE])
  raster_values <- popmaps_extract_tuning_values(prepared$rast, coords)
  distance_context <- popmaps_prepare_tuning_distance_context(
    surface = surface,
    surface_object = surface_object,
    prepared = prepared,
    coords = coords,
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
  axis_count <- ncol(locations) - 3
  validation_folds <- popmaps_make_validation_folds(
    locations = locations,
    validation = validation,
    n_blocks = n_blocks,
    block_assignments = block_assignments,
    spatial_block_repeats = spatial_block_repeats,
    spatial_block_seed = spatial_block_seed
  )

  scored_sites_per_grid <- sum(vapply(validation_folds, function(fold) {
    length(fold$assessment_idx)
  }, integer(1)))
  fold_rows <- vector("list", nrow(parameter_grid) * scored_sites_per_grid)
  fold_idx <- 1

  for (combo_idx in seq_len(nrow(parameter_grid))) {
    combo <- parameter_grid[combo_idx, ]

    for (validation_fold in validation_folds) {
      for (site_idx in validation_fold$assessment_idx) {
        prediction <- popmaps_predict_site_tuning(
          site_idx = site_idx,
          training_idx = validation_fold$analysis_idx,
          locations = locations,
          raster_value = raster_values[site_idx],
          threshold = threshold,
          empirical_pt_dist = combo$empirical_pt_dist,
          num_sites = combo$num_sites,
          num_tested = combo$num_tested,
          popmod = combo$popmod,
          dist_prob_func = dist_prob_func,
          distance_context = distance_context
        )

        fold_rows[[fold_idx]] <- popmaps_tuning_fold_row(
          combo = combo,
          validation = validation,
          validation_fold = validation_fold,
          site_idx = site_idx,
          site = as.character(locations$V1[site_idx]),
          prediction = prediction,
          axis_count = axis_count,
          distance_units = distance_context$distance_units
        )
        fold_idx <- fold_idx + 1
      }
    }
  }

  folds <- do.call(rbind, fold_rows)
  rownames(folds) <- NULL

  results <- popmaps_summarize_tuning_results(folds)
  best <- popmaps_tuning_best_row(results, primary_metric)

  tuning <- list(
    results = results,
    folds = folds,
    best = best,
    primary_metric = primary_metric,
    validation = validation,
    call = call
  )
  if (!is.null(search)) {
    tuning$search <- search
  }
  class(tuning) <- class

  if (!isTRUE(quiet)) {
    message(
      "Evaluated ",
      nrow(results),
      " parameter combinations across ",
      length(validation_folds),
      " validation folds."
    )
  }

  tuning
}

#' @export
print.popmaps_tuning <- function(x, ...) {
  cat("POPMAPS parameter tuning\n")
  cat("Primary metric: ", x$primary_metric, "\n", sep = "")
  cat("Validation: ", x$validation, "\n", sep = "")
  cat("Combinations: ", nrow(x$results), "\n", sep = "")
  cat("Folds: ", nrow(x$folds), "\n\n", sep = "")
  cat("Best parameter combination:\n")
  print(x$best, row.names = FALSE)
  invisible(x)
}

#' @export
print.popmaps_adaptive_tuning <- function(x, ...) {
  cat("Adaptive POPMAPS parameter tuning\n")
  cat("Method: ", x$search$method, "\n", sep = "")
  cat("Primary metric: ", x$primary_metric, "\n", sep = "")
  cat("Validation: ", x$validation, "\n", sep = "")
  cat("Combinations: ", nrow(x$results), "\n", sep = "")
  cat("Folds: ", nrow(x$folds), "\n\n", sep = "")
  cat("Best parameter combination:\n")
  print(x$best, row.names = FALSE)
  invisible(x)
}

#' @export
print.popmaps_tuning_grid <- function(x, ...) {
  distance_units <- if (!is.null(x$distance_units)) x$distance_units else "km"
  cat("Suggested POPMAPS tuning grid\n")
  cat("Distance units: ", distance_units, "\n", sep = "")
  cat("Reference distance: ", round(x$reference_distance, 3), " (", x$distance_reference, ")\n", sep = "")
  cat("num_sites: ", paste(x$num_sites, collapse = ", "), "\n", sep = "")
  cat("num_tested: ", paste(x$num_tested, collapse = ", "), "\n", sep = "")
  cat("empirical_pt_dist: ", paste(round(x$empirical_pt_dist, 3), collapse = ", "), "\n", sep = "")
  cat("popmod: ", paste(signif(x$popmod, 4), collapse = ", "), "\n", sep = "")
  invisible(x)
}

popmaps_default_num_sites <- function(n_training) {
  popmaps_check_whole_count(n_training, "`n_training`", positive = TRUE)

  candidates <- unique(round(c(0.33, 0.50, 0.75, 1.00) * n_training))
  candidates <- sort(unique(c(min(5, n_training), candidates)))
  candidates[candidates >= 2 & candidates <= n_training]
}

popmaps_make_validation_folds <- function(locations,
                                          validation,
                                          n_blocks,
                                          block_assignments = NULL,
                                          spatial_block_repeats = 1,
                                          spatial_block_seed = NULL) {
  n_sites <- nrow(locations)
  popmaps_check_whole_count(spatial_block_repeats, "`spatial_block_repeats`", positive = TRUE)
  if (!is.null(spatial_block_seed)) {
    popmaps_check_whole_count(spatial_block_seed, "`spatial_block_seed`", allow_zero = TRUE)
  }

  if (validation == "loo") {
    return(lapply(seq_len(n_sites), function(site_idx) {
      list(
        repeat_id = 1L,
        fold_id = site_idx,
        block_id = NA_character_,
        assessment_idx = site_idx,
        analysis_idx = setdiff(seq_len(n_sites), site_idx)
      )
    }))
  }

  if (validation != "spatial_block") {
    stop("Unsupported validation mode.", call. = FALSE)
  }

  if (!is.null(block_assignments) && spatial_block_repeats > 1) {
    stop("Repeated spatial-block validation requires `block_assignments = NULL`.", call. = FALSE)
  }

  old_seed <- NULL
  if (is.null(block_assignments) && spatial_block_repeats > 1 && !is.null(spatial_block_seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(spatial_block_seed)
  }

  if (is.null(block_assignments)) {
    popmaps_check_whole_count(n_blocks, "`n_blocks`", positive = TRUE)
    if (n_blocks < 2) {
      stop("`n_blocks` must be at least 2 for spatial-block validation.", call. = FALSE)
    }
    if (n_blocks > n_sites) {
      stop("`n_blocks` cannot exceed the number of empirical sites.", call. = FALSE)
    }
  } else {
    if (length(block_assignments) != n_sites) {
      stop("`block_assignments` must have one value per empirical site.", call. = FALSE)
    }
    if (any(is.na(block_assignments))) {
      stop("`block_assignments` cannot contain missing values.", call. = FALSE)
    }
    block_assignments <- as.character(block_assignments)
  }

  validation_folds <- list()
  fold_idx <- 1L
  for (repeat_id in seq_len(spatial_block_repeats)) {
    repeat_assignments <- if (is.null(block_assignments)) {
      angle <- if (repeat_id == 1L) 0 else stats::runif(1, min = 0, max = pi)
      popmaps_assign_spatial_blocks(locations, n_blocks, angle = angle)
    } else {
      block_assignments
    }

    block_ids <- unique(repeat_assignments)
    if (length(block_ids) < 2) {
      stop("Spatial-block validation requires at least two non-empty blocks.", call. = FALSE)
    }

    for (block_id in block_ids) {
      assessment_idx <- which(repeat_assignments == block_id)
      fold_block_id <- if (spatial_block_repeats == 1L) {
        block_id
      } else {
        paste0("repeat", repeat_id, "_", block_id)
      }
      validation_folds[[fold_idx]] <- list(
        repeat_id = repeat_id,
        fold_id = fold_idx,
        block_id = fold_block_id,
        assessment_idx = assessment_idx,
        analysis_idx = setdiff(seq_len(n_sites), assessment_idx)
      )
      fold_idx <- fold_idx + 1L
    }
  }

  validation_folds
}

popmaps_assign_spatial_blocks <- function(locations, n_blocks, angle = 0) {
  n_x <- ceiling(sqrt(n_blocks))
  n_y <- ceiling(n_blocks / n_x)
  coords <- as.matrix(locations[, c("V2", "V3"), drop = FALSE])
  coords <- scale(coords, center = TRUE, scale = TRUE)
  coords[!is.finite(coords)] <- 0
  rotated_x <- coords[, 1] * cos(angle) - coords[, 2] * sin(angle)
  rotated_y <- coords[, 1] * sin(angle) + coords[, 2] * cos(angle)

  x_bin <- popmaps_rank_bins(rotated_x, n_x)
  y_bin <- popmaps_rank_bins(rotated_y, n_y)
  paste0("x", x_bin, "_y", y_bin)
}

popmaps_rank_bins <- function(values, n_bins) {
  if (n_bins <= 1) {
    return(rep(1L, length(values)))
  }

  ranks <- rank(values, ties.method = "first")
  bins <- ceiling(ranks / length(values) * n_bins)
  pmax(1L, pmin(n_bins, bins))
}

popmaps_check_whole_count <- function(x,
                                      label,
                                      positive = FALSE,
                                      allow_zero = FALSE) {
  popmaps_check_finite_scalar(x, label)
  if (x != floor(x)) {
    stop(label, " must be a whole number.", call. = FALSE)
  }
  if (isTRUE(positive) && x < 1) {
    stop(label, " must be positive.", call. = FALSE)
  }
  if (isTRUE(allow_zero) && x < 0) {
    stop(label, " must be non-negative.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_make_tuning_grid <- function(num_sites,
                                     num_tested,
                                     popmod,
                                     empirical_pt_dist) {
  parameter_grid <- expand.grid(
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    empirical_pt_dist = empirical_pt_dist,
    KEEP.OUT.ATTRS = FALSE
  )

  popmaps_normalize_tuning_grid(parameter_grid)
}

popmaps_normalize_tuning_grid <- function(parameter_grid) {
  required_cols <- c("num_sites", "num_tested", "popmod", "empirical_pt_dist")
  if (!is.data.frame(parameter_grid) || !all(required_cols %in% names(parameter_grid))) {
    stop(
      "`parameter_grid` must contain num_sites, num_tested, popmod, and empirical_pt_dist columns.",
      call. = FALSE
    )
  }

  parameter_grid <- unique(parameter_grid[, required_cols])
  parameter_grid$num_sites <- as.integer(parameter_grid$num_sites)
  parameter_grid$num_tested <- as.integer(parameter_grid$num_tested)
  parameter_grid$popmod <- as.numeric(parameter_grid$popmod)
  parameter_grid$empirical_pt_dist <- as.numeric(parameter_grid$empirical_pt_dist)

  if (nrow(parameter_grid) < 1) {
    stop("At least one parameter combination is required.", call. = FALSE)
  }
  if (
    any(!is.finite(parameter_grid$num_sites)) ||
      any(parameter_grid$num_sites < 1) ||
      any(parameter_grid$num_sites != floor(parameter_grid$num_sites))
  ) {
    stop("`num_sites` must contain positive whole numbers.", call. = FALSE)
  }
  if (
    any(!is.finite(parameter_grid$num_tested)) ||
      any(parameter_grid$num_tested < 1) ||
      any(parameter_grid$num_tested != floor(parameter_grid$num_tested))
  ) {
    stop("`num_tested` must contain positive whole numbers.", call. = FALSE)
  }
  if (any(!is.finite(parameter_grid$popmod))) {
    stop("`popmod` must contain finite numeric values.", call. = FALSE)
  }
  if (
    any(!is.finite(parameter_grid$empirical_pt_dist)) ||
      any(parameter_grid$empirical_pt_dist < 0)
  ) {
    stop("`empirical_pt_dist` must contain non-negative finite numeric values.", call. = FALSE)
  }
  if (any(parameter_grid$num_tested > parameter_grid$num_sites)) {
    stop("Every tested parameter combination must have `num_tested <= num_sites`.", call. = FALSE)
  }

  parameter_grid$combo_id <- seq_len(nrow(parameter_grid))
  parameter_grid[, c("combo_id", required_cols)]
}

popmaps_resolve_search_values <- function(values,
                                          default_values,
                                          label,
                                          whole_number = FALSE,
                                          positive = FALSE,
                                          nonnegative = FALSE) {
  if (is.null(values)) {
    values <- default_values
  }

  popmaps_check_tuning_values(
    values,
    label,
    whole_number = whole_number,
    positive = positive,
    nonnegative = nonnegative
  )
}

popmaps_sample_tuning_space <- function(n, search_space, method) {
  popmaps_check_whole_count(n, "`n`", positive = TRUE)

  if (min(search_space$num_tested) > max(search_space$num_sites)) {
    stop("The search space has no possible combinations with `num_tested <= num_sites`.", call. = FALSE)
  }

  parameter_grid <- data.frame()
  attempts <- 0
  while (nrow(parameter_grid) < n && attempts < 20) {
    attempts <- attempts + 1
    batch_n <- max(n - nrow(parameter_grid), n)
    sampled_sites <- popmaps_sample_discrete_values(search_space$num_sites, batch_n, method)
    sampled_tested <- vapply(sampled_sites, function(site_count) {
      valid_tested <- search_space$num_tested[search_space$num_tested <= site_count]
      if (length(valid_tested) < 1) {
        return(NA_real_)
      }
      sample(valid_tested, 1)
    }, numeric(1))

    batch <- data.frame(
      num_sites = sampled_sites,
      num_tested = sampled_tested,
      popmod = popmaps_sample_numeric_values(search_space$popmod, batch_n, method),
      empirical_pt_dist = popmaps_sample_numeric_values(
        search_space$empirical_pt_dist,
        batch_n,
        method
      )
    )
    batch <- batch[is.finite(batch$num_tested), , drop = FALSE]
    parameter_grid <- unique(rbind(parameter_grid, batch))
  }

  if (nrow(parameter_grid) < 1) {
    stop("No valid parameter combinations could be sampled.", call. = FALSE)
  }
  if (nrow(parameter_grid) > n) {
    parameter_grid <- parameter_grid[seq_len(n), , drop = FALSE]
  }

  popmaps_normalize_tuning_grid(parameter_grid)
}

popmaps_sample_discrete_values <- function(values, n, method) {
  values <- sort(unique(values))
  if (length(values) == 1) {
    return(rep(values, n))
  }
  if (method == "random") {
    return(sample(values, n, replace = TRUE))
  }

  u <- sample((seq_len(n) - stats::runif(n)) / n)
  idx <- pmax(1, pmin(length(values), ceiling(u * length(values))))
  values[idx]
}

popmaps_sample_numeric_values <- function(values, n, method) {
  values <- sort(unique(values))
  if (length(values) == 1) {
    return(rep(values, n))
  }

  value_range <- range(values)
  u <- if (method == "random") {
    stats::runif(n)
  } else {
    sample((seq_len(n) - stats::runif(n)) / n)
  }
  value_range[1] + (u * diff(value_range))
}

popmaps_refine_tuning_space <- function(results, search_space, primary_metric) {
  ordered_results <- results[popmaps_order_tuning_results(results, primary_metric), , drop = FALSE]
  top_n <- max(3, ceiling(nrow(ordered_results) * 0.15))
  top_results <- ordered_results[seq_len(min(top_n, nrow(ordered_results))), , drop = FALSE]

  list(
    num_sites = popmaps_neighbor_values(search_space$num_sites, top_results$num_sites),
    num_tested = popmaps_neighbor_values(search_space$num_tested, top_results$num_tested),
    popmod = popmaps_refine_numeric_values(top_results$popmod, search_space$popmod),
    empirical_pt_dist = popmaps_refine_numeric_values(
      top_results$empirical_pt_dist,
      search_space$empirical_pt_dist,
      lower_bound = 0
    )
  )
}

popmaps_neighbor_values <- function(values, centers) {
  values <- sort(unique(values))
  center_idx <- match(unique(centers), values)
  center_idx <- center_idx[!is.na(center_idx)]
  neighbor_idx <- unique(unlist(lapply(center_idx, function(idx) seq.int(idx - 1, idx + 1))))
  neighbor_idx <- neighbor_idx[neighbor_idx >= 1 & neighbor_idx <= length(values)]
  values[neighbor_idx]
}

popmaps_refine_numeric_values <- function(top_values, full_values, lower_bound = -Inf) {
  full_range <- range(full_values)
  if (length(unique(full_values)) == 1 || diff(full_range) == 0) {
    return(unique(full_values))
  }

  top_range <- range(top_values)
  top_span <- diff(top_range)
  if (top_span == 0) {
    top_span <- diff(full_range) * 0.20
  }

  refined <- c(top_range[1] - top_span, top_range[2] + top_span)
  refined[1] <- max(refined[1], full_range[1], lower_bound)
  refined[2] <- min(refined[2], full_range[2])
  if (refined[1] == refined[2]) {
    refined <- refined[1]
  }

  refined
}

popmaps_metric_is_maximized <- function(primary_metric) {
  primary_metric %in% c("dominant_accuracy", "dominant_probability")
}

popmaps_percent_change <- function(delta, reference) {
  if (!is.finite(delta) || !is.finite(reference) || abs(reference) < .Machine$double.eps) {
    return(NA_real_)
  }

  100 * delta / abs(reference)
}

popmaps_order_tuning_results <- function(results, primary_metric) {
  maximize <- popmaps_metric_is_maximized(primary_metric)
  metric_values <- results[[primary_metric]]
  order_values <- if (maximize) -metric_values else metric_values
  order_values[!is.finite(order_values)] <- Inf
  metric_sd_col <- paste0(primary_metric, "_repeat_sd")
  metric_sd_values <- if (metric_sd_col %in% names(results)) {
    results[[metric_sd_col]]
  } else {
    rep(NA_real_, nrow(results))
  }
  metric_sd_values[!is.finite(metric_sd_values)] <- Inf

  order(
    results$failed_folds,
    -results$n_scored,
    order_values,
    metric_sd_values,
    results$num_tested,
    results$num_sites,
    results$popmod,
    results$empirical_pt_dist
  )
}

popmaps_decay_distance <- function(popmod, retained_weight) {
  distance <- rep(NA_real_, length(popmod))
  no_decay <- popmod == 0
  distance[no_decay] <- Inf
  decay <- popmod < 0
  distance[decay] <- log(retained_weight) / popmod[decay]
  distance
}

popmaps_check_tuning_values <- function(x,
                                        label,
                                        whole_number = FALSE,
                                        positive = FALSE,
                                        nonnegative = FALSE) {
  if (!is.numeric(x) || length(x) < 1 || any(!is.finite(x))) {
    stop(label, " must contain finite numeric values.", call. = FALSE)
  }
  if (isTRUE(whole_number) && any(x != floor(x))) {
    stop(label, " must contain whole numbers.", call. = FALSE)
  }
  if (isTRUE(positive) && any(x < 1)) {
    stop(label, " must contain positive values.", call. = FALSE)
  }
  if (isTRUE(nonnegative) && any(x < 0)) {
    stop(label, " must contain non-negative values.", call. = FALSE)
  }

  sort(unique(x))
}

popmaps_extract_tuning_values <- function(rast, coords) {
  extracted <- terra::extract(rast, coords)
  if (is.data.frame(extracted) || is.matrix(extracted)) {
    return(as.numeric(extracted[, ncol(extracted)]))
  }

  as.numeric(extracted)
}

popmaps_prepare_tuning_distance_context <- function(surface,
                                                    surface_object,
                                                    prepared,
                                                    coords,
                                                    surface_values,
                                                    rescale_conductance,
                                                    resistance_epsilon) {
  if (surface == "G") {
    return(list(
      surface = "G",
      distance_units = "km",
      site_distances = popmaps_empirical_site_distances(coords)
    ))
  }

  if (is.null(surface_object)) {
    surface_object <- prepare_popmaps_surface(
      input_raster = prepared$rast,
      surface = "C",
      surface_values = surface_values,
      rescale_conductance = rescale_conductance,
      resistance_epsilon = resistance_epsilon
    )
  }

  graph <- popmaps_cost_distance_graph(surface_object, directions = 8)
  site_distances <- tryCatch(
    popmaps_cost_distance_matrix(
      surface = surface_object,
      from_coords = coords,
      directions = 8,
      graph = graph
    ),
    error = function(err) {
      stop(
        "Could not calculate least-cost distances for empirical sites: ",
        conditionMessage(err),
        call. = FALSE
      )
    }
  )

  list(
    surface = "C",
    distance_units = "cost_distance",
    site_distances = site_distances,
    graph = graph,
    surface_object = surface_object
  )
}

popmaps_predict_site_tuning <- function(site_idx,
                                        training_idx = NULL,
                                        locations,
                                        raster_value,
                                        threshold,
                                        empirical_pt_dist,
                                        num_sites,
                                        num_tested,
                                        popmod,
                                        dist_prob_func,
                                        distance_context) {
  if (is.null(training_idx)) {
    training_idx <- setdiff(seq_len(nrow(locations)), site_idx)
  }

  site_distances <- distance_context$site_distances[site_idx, training_idx]
  empirical_distances <- distance_context$site_distances[training_idx, training_idx, drop = FALSE]

  popmaps_predict_site_from_distances(
    site_idx = site_idx,
    training_idx = training_idx,
    locations = locations,
    raster_value = raster_value,
    threshold = threshold,
    empirical_pt_dist = empirical_pt_dist,
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    dist_prob_func = dist_prob_func,
    site_distances = site_distances,
    empirical_distances = empirical_distances
  )
}

popmaps_predict_site_geographic <- function(site_idx,
                                            training_idx = NULL,
                                            locations,
                                            raster_value,
                                            threshold,
                                            empirical_pt_dist,
                                            num_sites,
                                            num_tested,
                                            popmod,
                                            dist_prob_func) {
  if (is.null(training_idx)) {
    training_idx <- setdiff(seq_len(nrow(locations)), site_idx)
  }
  training <- locations[training_idx, , drop = FALSE]
  training_coords <- as.matrix(training[, 2:3, drop = FALSE])
  site_distances <- popmaps_earth_dist(
    lat1 = locations$V3[site_idx],
    long1 = locations$V2[site_idx],
    lat2 = training$V3,
    long2 = training$V2
  )
  empirical_distances <- popmaps_empirical_site_distances(training_coords)

  popmaps_predict_site_from_distances(
    site_idx = site_idx,
    training_idx = training_idx,
    locations = locations,
    raster_value = raster_value,
    threshold = threshold,
    empirical_pt_dist = empirical_pt_dist,
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    dist_prob_func = dist_prob_func,
    site_distances = site_distances,
    empirical_distances = empirical_distances
  )
}

popmaps_predict_site_from_distances <- function(site_idx,
                                                training_idx,
                                                locations,
                                                raster_value,
                                                threshold,
                                                empirical_pt_dist,
                                                num_sites,
                                                num_tested,
                                                popmod,
                                                dist_prob_func,
                                                site_distances,
                                                empirical_distances) {
  axis_cols <- seq.int(4, ncol(locations))
  observed <- popmaps_normalize_probability(as.numeric(locations[site_idx, axis_cols]))

  if (any(is.na(observed))) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "Observed ancestry coefficients could not be normalized."
    ))
  }
  if (is.na(raster_value)) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "Withheld site falls outside non-NA raster cells."
    ))
  }
  if (raster_value < threshold) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "Withheld site raster value is below `threshold`."
    ))
  }
  if (length(training_idx) < num_tested) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "Validation fold has fewer training sites than `num_tested`."
    ))
  }
  if (length(training_idx) < num_sites) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "Validation fold has fewer training sites than `num_sites`."
    ))
  }

  reachable_sites <- which(is.finite(site_distances))
  if (length(reachable_sites) < num_sites) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "Validation fold has fewer reachable training sites than `num_sites`."
    ))
  }

  training <- locations[training_idx, , drop = FALSE]
  training_ancestry <- as.matrix(training[, axis_cols, drop = FALSE])

  candidate_sites <- reachable_sites[order(site_distances[reachable_sites])][seq_len(num_sites)]

  selected_sites <- tryCatch(
    popmaps_select_empirical_sites(
      candidate_sites = candidate_sites,
      empirical_distances = empirical_distances,
      empirical_pt_dist = empirical_pt_dist,
      num_tested = num_tested
    ),
    error = function(err) {
      err
    }
  )
  if (inherits(selected_sites, "error")) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = conditionMessage(selected_sites)
    ))
  }

  site_weights <- vapply(
    site_distances[selected_sites],
    function(distance) dist_prob_func(popmod, distance),
    numeric(1)
  )

  if (any(!is.finite(site_weights))) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "`dist_prob_func` returned non-finite values."
    ))
  }

  cell_prob <- colSums(training_ancestry[selected_sites, , drop = FALSE] * (site_weights / num_tested))
  predicted <- popmaps_normalize_probability(cell_prob)
  if (any(is.na(predicted))) {
    return(popmaps_failed_tuning_prediction(
      observed = observed,
      message = "Predicted ancestry coefficients could not be normalized."
    ))
  }

  list(
    predicted = predicted,
    observed = observed,
    metrics = popmaps_score_prediction(predicted, observed),
    message = NA_character_
  )
}

popmaps_failed_tuning_prediction <- function(observed, message) {
  list(
    predicted = rep(NA_real_, length(observed)),
    observed = observed,
    metrics = c(
      mae = NA_real_,
      rmse = NA_real_,
      hellinger = NA_real_,
      dominant_accuracy = NA_real_,
      dominant_probability = NA_real_
    ),
    message = message
  )
}

popmaps_normalize_probability <- function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < 0)) {
    return(rep(NA_real_, length(x)))
  }

  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(x)))
  }

  x / total
}

popmaps_score_prediction <- function(predicted, observed) {
  diff <- predicted - observed

  c(
    mae = mean(abs(diff)),
    rmse = sqrt(mean(diff^2)),
    hellinger = sqrt(sum((sqrt(predicted) - sqrt(observed))^2)) / sqrt(2),
    dominant_accuracy = as.numeric(which.max(predicted) == which.max(observed)),
    dominant_probability = predicted[which.max(observed)]
  )
}

popmaps_tuning_fold_row <- function(combo,
                                    validation,
                                    validation_fold,
                                    site_idx,
                                    site,
                                    prediction,
                                    axis_count,
                                    distance_units = "km") {
  predicted_names <- paste0("predicted_axis_", seq_len(axis_count))
  observed_names <- paste0("observed_axis_", seq_len(axis_count))

  cbind(
    data.frame(
      combo_id = combo$combo_id,
      validation = validation,
      repeat_id = validation_fold$repeat_id,
      fold_id = validation_fold$fold_id,
      block_id = validation_fold$block_id,
      site_index = site_idx,
      site = site,
      n_training = length(validation_fold$analysis_idx),
      num_sites = combo$num_sites,
      num_tested = combo$num_tested,
      popmod = combo$popmod,
      half_distance = popmaps_decay_distance(combo$popmod, 0.5),
      ten_pct_distance = popmaps_decay_distance(combo$popmod, 0.1),
      half_distance_km = popmaps_decay_distance(combo$popmod, 0.5),
      ten_pct_distance_km = popmaps_decay_distance(combo$popmod, 0.1),
      distance_units = distance_units,
      empirical_pt_dist = combo$empirical_pt_dist,
      mae = prediction$metrics[["mae"]],
      rmse = prediction$metrics[["rmse"]],
      hellinger = prediction$metrics[["hellinger"]],
      dominant_accuracy = prediction$metrics[["dominant_accuracy"]],
      dominant_probability = prediction$metrics[["dominant_probability"]],
      message = prediction$message,
      stringsAsFactors = FALSE
    ),
    stats::setNames(
      as.data.frame(as.list(prediction$predicted), stringsAsFactors = FALSE),
      predicted_names
    ),
    stats::setNames(
      as.data.frame(as.list(prediction$observed), stringsAsFactors = FALSE),
      observed_names
    )
  )
}

popmaps_summarize_tuning_results <- function(folds) {
  metric_cols <- c(
    "mae",
    "rmse",
    "hellinger",
    "dominant_accuracy",
    "dominant_probability"
  )
  combo_ids <- unique(folds$combo_id)

  result_rows <- lapply(combo_ids, function(combo_id) {
    idx <- folds$combo_id == combo_id
    first <- folds[idx, ][1, ]

    metrics <- vapply(metric_cols, function(metric) {
      value <- mean(folds[[metric]][idx], na.rm = TRUE)
      if (is.nan(value)) {
        NA_real_
      } else {
        value
      }
    }, numeric(1))
    repeat_sd <- vapply(metric_cols, function(metric) {
      repeat_means <- tapply(
        folds[[metric]][idx],
        folds$repeat_id[idx],
        mean,
        na.rm = TRUE
      )
      repeat_means[is.nan(repeat_means)] <- NA_real_
      repeat_means <- repeat_means[is.finite(repeat_means)]
      if (length(repeat_means) > 1) {
        stats::sd(repeat_means)
      } else {
        NA_real_
      }
    }, numeric(1))
    names(repeat_sd) <- paste0(metric_cols, "_repeat_sd")

    data.frame(
      combo_id = combo_id,
      validation = first$validation,
      num_sites = first$num_sites,
      num_tested = first$num_tested,
      popmod = first$popmod,
      half_distance = first$half_distance,
      ten_pct_distance = first$ten_pct_distance,
      half_distance_km = first$half_distance_km,
      ten_pct_distance_km = first$ten_pct_distance_km,
      distance_units = first$distance_units,
      empirical_pt_dist = first$empirical_pt_dist,
      n_validation_repeats = length(unique(folds$repeat_id[idx])),
      n_validation_folds = length(unique(folds$fold_id[idx])),
      n_folds = sum(idx),
      n_scored = sum(is.finite(folds$rmse[idx])),
      n_training_min = min(folds$n_training[idx]),
      failed_folds = sum(!is.na(folds$message[idx]) & nzchar(folds$message[idx])),
      t(metrics),
      t(repeat_sd),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, result_rows)
  rownames(results) <- NULL
  results
}

popmaps_tuning_best_row <- function(results, primary_metric) {
  best_idx <- popmaps_order_tuning_results(results, primary_metric)[1]

  results[best_idx, , drop = FALSE]
}
