#!/usr/bin/env Rscript

truthy_env <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

read_arg_or_env <- function(args, index, env, required = TRUE, default = NULL) {
  value <- if (length(args) >= index && nzchar(args[[index]])) {
    args[[index]]
  } else {
    Sys.getenv(env, unset = if (is.null(default)) "" else default)
  }

  if (required && (is.null(value) || is.na(value) || !nzchar(value))) {
    stop("Missing value. Supply command argument ", index, " or set ", env, ".", call. = FALSE)
  }

  value
}

read_integer_env <- function(env, default, allow_zero = FALSE) {
  value <- as.numeric(Sys.getenv(env, unset = as.character(default)))
  if (
    length(value) != 1 ||
      is.na(value) ||
      !is.finite(value) ||
      value != floor(value) ||
      value < if (allow_zero) 0 else 1
  ) {
    stop(env, " must be a ", if (allow_zero) "non-negative" else "positive", " whole number.", call. = FALSE)
  }

  as.integer(value)
}

read_modes_env <- function(env, default) {
  modes <- strsplit(Sys.getenv(env, unset = default), ",", fixed = TRUE)[[1]]
  modes <- trimws(modes)
  modes <- modes[nzchar(modes)]
  allowed <- c("loo", "spatial_block")
  if (length(modes) < 1 || any(!modes %in% allowed)) {
    stop(env, " must contain one or more of: ", paste(allowed, collapse = ", "), call. = FALSE)
  }

  unique(modes)
}

find_example_pairs <- function(input_dir) {
  rasters <- list.files(input_dir, pattern = "_avg[.]asc$", full.names = TRUE)
  rows <- lapply(rasters, function(raster_path) {
    species <- sub("_avg[.]asc$", "", basename(raster_path))
    locations_path <- file.path(dirname(raster_path), paste0(species, ".txt"))
    if (!file.exists(locations_path)) {
      return(NULL)
    }

    data.frame(
      species = species,
      raster_path = raster_path,
      locations_path = locations_path,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) < 1) {
    stop("No '*_avg.asc' and matching '*.txt' example pairs found in ", input_dir, ".", call. = FALSE)
  }

  do.call(rbind, rows)
}

load_popmaps2 <- function() {
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("popmaps2", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  stop(
    "The popmaps2 package must be installed before running this script. ",
    "From the repository root, run: R CMD INSTALL .",
    call. = FALSE
  )
}

write_table <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  message("Wrote ", path)
}

args <- commandArgs(trailingOnly = TRUE)
default_input_dir <- file.path(dirname(getwd()), "popmaps_test_data")
input_dir <- read_arg_or_env(args, 1, "POPMAPS_EXAMPLE_DIR", default = default_input_dir)
output_dir <- read_arg_or_env(
  args,
  2,
  "POPMAPS_EXAMPLE_OUTPUT_DIR",
  required = FALSE,
  default = file.path(input_dir, "tuning_outputs")
)

load_popmaps2()
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("The terra package is required for empirical tuning validation.", call. = FALSE)
}

input_dir <- normalizePath(input_dir, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

search <- Sys.getenv("POPMAPS_EXAMPLE_SEARCH", unset = "grid")
if (!search %in% c("grid", "adaptive")) {
  stop("POPMAPS_EXAMPLE_SEARCH must be 'grid' or 'adaptive'.", call. = FALSE)
}
validation_modes <- read_modes_env("POPMAPS_EXAMPLE_VALIDATION", "loo,spatial_block")
aggregate_fact <- read_integer_env("POPMAPS_EXAMPLE_AGGREGATE", 1)
n_blocks <- read_integer_env("POPMAPS_EXAMPLE_N_BLOCKS", 4)
n_initial <- read_integer_env("POPMAPS_EXAMPLE_N_INITIAL", 50)
n_refine <- read_integer_env("POPMAPS_EXAMPLE_N_REFINE", 50, allow_zero = TRUE)
near_best_tolerance <- as.numeric(Sys.getenv("POPMAPS_EXAMPLE_NEAR_BEST_TOLERANCE", unset = "0.05"))
if (length(near_best_tolerance) != 1 || !is.finite(near_best_tolerance) || near_best_tolerance < 0) {
  stop("POPMAPS_EXAMPLE_NEAR_BEST_TOLERANCE must be one non-negative number.", call. = FALSE)
}
seed <- read_integer_env("POPMAPS_EXAMPLE_SEED", 1, allow_zero = TRUE)
write_folds <- truthy_env(Sys.getenv("POPMAPS_EXAMPLE_WRITE_FOLDS", unset = "false"))

pairs <- find_example_pairs(input_dir)
message("Found ", nrow(pairs), " empirical example datasets.")

run_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
summary_rows <- list()
overview_rows <- list()
range_rows <- list()
effect_rows <- list()

for (row_idx in seq_len(nrow(pairs))) {
  species <- pairs$species[row_idx]
  message("Reading ", species)
  raster_surface <- terra::rast(pairs$raster_path[row_idx])
  if (aggregate_fact > 1) {
    raster_surface <- terra::aggregate(raster_surface, fact = aggregate_fact, fun = mean, na.rm = TRUE)
  }
  locations <- utils::read.table(
    pairs$locations_path[row_idx],
    header = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- popmaps2::suggest_tuning_grid(locations)

  for (validation in validation_modes) {
    message("Tuning ", species, " with validation = ", validation, " and search = ", search)
    tuning <- if (search == "adaptive") {
      popmaps2::adaptive_tune_popmaps(
        input_raster = raster_surface,
        input_locs = locations,
        validation = validation,
        n_blocks = n_blocks,
        n_initial = n_initial,
        n_refine = n_refine,
        seed = seed,
        quiet = TRUE
      )
    } else {
      popmaps2::tune_popmaps(
        input_raster = raster_surface,
        input_locs = locations,
        validation = validation,
        n_blocks = n_blocks,
        empirical_pt_dist = grid$empirical_pt_dist,
        num_sites = grid$num_sites,
        num_tested = grid$num_tested,
        popmod = grid$popmod,
        quiet = TRUE
      )
    }

    diagnostics <- popmaps2::diagnose_tuning(
      tuning,
      near_best_tolerance = near_best_tolerance
    )
    prefix <- paste(species, validation, search, run_stamp, sep = "-")

    write_table(tuning$results, file.path(output_dir, paste0(prefix, "-results.csv")))
    if (write_folds) {
      write_table(tuning$folds, file.path(output_dir, paste0(prefix, "-folds.csv")))
    }
    write_table(diagnostics$near_best, file.path(output_dir, paste0(prefix, "-near-best.csv")))
    write_table(diagnostics$parameter_ranges, file.path(output_dir, paste0(prefix, "-parameter-ranges.csv")))
    write_table(diagnostics$parameter_effects, file.path(output_dir, paste0(prefix, "-parameter-effects.csv")))

    overview <- cbind(
      data.frame(species = species, search = search, stringsAsFactors = FALSE),
      diagnostics$overview
    )
    best <- cbind(
      data.frame(species = species, search = search, stringsAsFactors = FALSE),
      tuning$best
    )
    overview_rows[[length(overview_rows) + 1]] <- overview
    summary_rows[[length(summary_rows) + 1]] <- best
    ranges <- cbind(data.frame(species = species, validation = validation, search = search), diagnostics$parameter_ranges)
    effects <- cbind(data.frame(species = species, validation = validation, search = search), diagnostics$parameter_effects)
    range_rows[[length(range_rows) + 1]] <- ranges
    effect_rows[[length(effect_rows) + 1]] <- effects
  }
}

summary_table <- do.call(rbind, summary_rows)
overview_table <- do.call(rbind, overview_rows)
range_table <- do.call(rbind, range_rows)
effect_table <- do.call(rbind, effect_rows)

write_table(summary_table, file.path(output_dir, paste0("empirical-tuning-best-", search, "-", run_stamp, ".csv")))
write_table(overview_table, file.path(output_dir, paste0("empirical-tuning-overview-", search, "-", run_stamp, ".csv")))
write_table(range_table, file.path(output_dir, paste0("empirical-tuning-parameter-ranges-", search, "-", run_stamp, ".csv")))
write_table(effect_table, file.path(output_dir, paste0("empirical-tuning-parameter-effects-", search, "-", run_stamp, ".csv")))

print(overview_table)
