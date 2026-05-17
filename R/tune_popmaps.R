#' Tune POPMAPS parameters with leave-one-site-out validation
#'
#' @description `tune_popmaps()` compares POPMAPS parameter combinations by
#'   withholding each empirical sampling site, predicting its ancestry
#'   coefficients from the remaining sites, and summarizing prediction error.
#'   This is a faster, more structured geographic-distance alternative to the
#'   legacy `jackknife()` workflow.
#'
#' @param input_raster A `terra::SpatRaster`, `raster::RasterLayer`, or path to
#'   a raster file defining the interpolation surface.
#' @param input_locs A data frame or matrix with sampling location name,
#'   longitude, latitude, and one or more ancestry coefficient columns.
#' @param surface Character. Currently only `"G"` is supported.
#' @param empirical_pt_dist Minimum geographic distance, in kilometers, required
#'   between empirical sites selected for a prediction.
#' @param num_sites Integer vector. Candidate pool sizes to evaluate.
#' @param num_tested Integer vector. Numbers of empirical sites used to estimate
#'   ancestry coefficients.
#' @param popmod Numeric vector. Distance-decay parameter values to evaluate.
#' @param threshold Numeric scalar. Raster values below this threshold are not
#'   scored.
#' @param primary_metric Metric used to select the best parameter combination.
#' @param dist_prob_func Function defining the relationship between distance and
#'   empirical-site contribution.
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
                         empirical_pt_dist = 5,
                         num_sites = 10,
                         num_tested = c(2, 3, 4, 5, 6, 7, 8),
                         popmod = c(-0.001, -0.01, -0.05, -0.1, -0.15),
                         threshold = 0,
                         primary_metric = c("rmse", "mae", "hellinger",
                                            "dominant_accuracy",
                                            "dominant_probability"),
                         dist_prob_func = function(popmod_temp, distance) {
                           exp(popmod_temp * distance)
                         },
                         quiet = TRUE) {
  surface <- match.arg(surface, c("G", "C"))
  primary_metric <- match.arg(primary_metric)

  if (surface != "G") {
    stop(
      "`tune_popmaps()` currently supports only geographic-distance tuning with surface = 'G'.",
      call. = FALSE
    )
  }
  if (!is.function(dist_prob_func)) {
    stop("`dist_prob_func` must be a function.", call. = FALSE)
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
  popmaps_check_finite_scalar(empirical_pt_dist, "`empirical_pt_dist`")
  if (empirical_pt_dist < 0) {
    stop("`empirical_pt_dist` must be non-negative.", call. = FALSE)
  }

  parameter_grid <- expand.grid(
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    KEEP.OUT.ATTRS = FALSE
  )
  parameter_grid$combo_id <- seq_len(nrow(parameter_grid))
  parameter_grid <- parameter_grid[, c("combo_id", "num_sites", "num_tested", "popmod")]

  if (any(parameter_grid$num_tested > parameter_grid$num_sites)) {
    stop("Every tested parameter combination must have `num_tested <= num_sites`.", call. = FALSE)
  }

  prepared <- popmaps_prepare_inputs(
    input_raster = input_raster,
    input_locs = input_locs,
    surface = surface,
    num_sites = max(num_sites),
    num_tested = max(num_tested),
    threshold = threshold,
    empirical_pt_dist = empirical_pt_dist,
    jackknife = TRUE
  )

  locations <- prepared$locations
  coords <- as.matrix(locations[, 2:3, drop = FALSE])
  raster_values <- popmaps_extract_tuning_values(prepared$rast, coords)
  axis_count <- ncol(locations) - 3

  fold_rows <- vector("list", nrow(parameter_grid) * nrow(locations))
  fold_idx <- 1

  for (combo_idx in seq_len(nrow(parameter_grid))) {
    combo <- parameter_grid[combo_idx, ]

    for (site_idx in seq_len(nrow(locations))) {
      prediction <- popmaps_predict_site_geographic(
        site_idx = site_idx,
        locations = locations,
        raster_value = raster_values[site_idx],
        threshold = threshold,
        empirical_pt_dist = empirical_pt_dist,
        num_sites = combo$num_sites,
        num_tested = combo$num_tested,
        popmod = combo$popmod,
        dist_prob_func = dist_prob_func
      )

      fold_rows[[fold_idx]] <- popmaps_tuning_fold_row(
        combo = combo,
        empirical_pt_dist = empirical_pt_dist,
        site_idx = site_idx,
        site = as.character(locations$V1[site_idx]),
        prediction = prediction,
        axis_count = axis_count
      )
      fold_idx <- fold_idx + 1
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
    call = match.call()
  )
  class(tuning) <- "popmaps_tuning"

  if (!isTRUE(quiet)) {
    message(
      "Evaluated ",
      nrow(results),
      " parameter combinations across ",
      nrow(locations),
      " leave-one-site-out folds."
    )
  }

  tuning
}

print.popmaps_tuning <- function(x, ...) {
  cat("POPMAPS parameter tuning\n")
  cat("Primary metric: ", x$primary_metric, "\n", sep = "")
  cat("Combinations: ", nrow(x$results), "\n", sep = "")
  cat("Folds: ", nrow(x$folds), "\n\n", sep = "")
  cat("Best parameter combination:\n")
  print(x$best, row.names = FALSE)
  invisible(x)
}

popmaps_check_tuning_values <- function(x,
                                        label,
                                        whole_number = FALSE,
                                        positive = FALSE) {
  if (!is.numeric(x) || length(x) < 1 || any(!is.finite(x))) {
    stop(label, " must contain finite numeric values.", call. = FALSE)
  }
  if (isTRUE(whole_number) && any(x != floor(x))) {
    stop(label, " must contain whole numbers.", call. = FALSE)
  }
  if (isTRUE(positive) && any(x < 1)) {
    stop(label, " must contain positive values.", call. = FALSE)
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

popmaps_predict_site_geographic <- function(site_idx,
                                            locations,
                                            raster_value,
                                            threshold,
                                            empirical_pt_dist,
                                            num_sites,
                                            num_tested,
                                            popmod,
                                            dist_prob_func) {
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

  training <- locations[-site_idx, , drop = FALSE]
  training_coords <- as.matrix(training[, 2:3, drop = FALSE])
  training_ancestry <- as.matrix(training[, axis_cols, drop = FALSE])

  site_distances <- popmaps_earth_dist(
    lat1 = locations$V3[site_idx],
    long1 = locations$V2[site_idx],
    lat2 = training$V3,
    long2 = training$V2
  )
  candidate_sites <- order(site_distances)[seq_len(num_sites)]
  empirical_distances <- popmaps_empirical_site_distances(training_coords)

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
                                    empirical_pt_dist,
                                    site_idx,
                                    site,
                                    prediction,
                                    axis_count) {
  predicted_names <- paste0("predicted_axis_", seq_len(axis_count))
  observed_names <- paste0("observed_axis_", seq_len(axis_count))

  cbind(
    data.frame(
      combo_id = combo$combo_id,
      site_index = site_idx,
      site = site,
      num_sites = combo$num_sites,
      num_tested = combo$num_tested,
      popmod = combo$popmod,
      empirical_pt_dist = empirical_pt_dist,
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

    data.frame(
      combo_id = combo_id,
      num_sites = first$num_sites,
      num_tested = first$num_tested,
      popmod = first$popmod,
      empirical_pt_dist = first$empirical_pt_dist,
      n_folds = sum(idx),
      n_scored = sum(is.finite(folds$rmse[idx])),
      failed_folds = sum(!is.na(folds$message[idx]) & nzchar(folds$message[idx])),
      t(metrics),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, result_rows)
  rownames(results) <- NULL
  results
}

popmaps_tuning_best_row <- function(results, primary_metric) {
  maximize <- primary_metric %in% c("dominant_accuracy", "dominant_probability")
  metric_values <- results[[primary_metric]]
  order_values <- if (maximize) -metric_values else metric_values
  order_values[!is.finite(order_values)] <- Inf

  best_idx <- order(
    order_values,
    results$failed_folds,
    results$num_tested,
    results$num_sites,
    results$popmod
  )[1]

  results[best_idx, , drop = FALSE]
}
