#!/usr/bin/env Rscript

# Local ASLO runtime and memory benchmark.
#
# Purpose:
# - run `popmaps()` across one or more raster aggregation factors;
# - compare `surface = "G"` and `surface = "C"` runtimes when requested;
# - record elapsed time, coarse R memory, raster size, and resource settings.
#
# Usage:
#   Rscript tools/benchmark-aslo-local.R aslo_avg.asc aslo.txt output_dir
#
# Key environment controls:
# - POPMAPS_BENCH_AGGREGATES: comma-separated aggregate factors, e.g. 16,4,1.
# - POPMAPS_BENCH_SURFACES: comma-separated surfaces, e.g. G,C.
# - POPMAPS_BENCH_NUM_SITES, POPMAPS_BENCH_NUM_TESTED, POPMAPS_BENCH_POPMOD:
#   interpolation parameters held constant across benchmark rows.

# Resolve the shared tools directory whether called from the repository root or
# another working directory.
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE))
} else {
  file.path(getwd(), "tools")
}
source(file.path(script_dir, "popmaps-script-utils.R"))
resource_config <- popmaps_configure_script_resources("POPMAPS_BENCH")

# Argument/environment readers fail early when benchmark settings are malformed;
# this avoids spending time on a run that cannot be compared later.
read_arg_or_env <- function(args, index, env, label, required = TRUE, default = NULL) {
  value <- if (length(args) >= index && nzchar(args[[index]])) {
    args[[index]]
  } else {
    Sys.getenv(env, unset = default)
  }

  if (required && (is.null(value) || is.na(value) || !nzchar(value))) {
    stop(
      "Missing ",
      label,
      ". Supply command argument ",
      index,
      " or set ",
      env,
      ".",
      call. = FALSE
    )
  }

  value
}

read_numeric_env <- function(env, default) {
  value <- as.numeric(Sys.getenv(env, unset = as.character(default)))
  if (length(value) != 1 || is.na(value) || !is.finite(value)) {
    stop(env, " must be one finite numeric value.", call. = FALSE)
  }

  value
}

read_integer_env <- function(env, default) {
  value <- read_numeric_env(env, default)
  if (value < 1 || value != floor(value)) {
    stop(env, " must be a positive whole number.", call. = FALSE)
  }

  as.integer(value)
}

# Read comma-separated controls such as aggregation factors and surface modes.
read_csv_values <- function(env, default, numeric = FALSE) {
  values <- strsplit(Sys.getenv(env, unset = default), ",", fixed = TRUE)[[1]]
  values <- trimws(values)
  values <- values[nzchar(values)]
  if (length(values) < 1) {
    stop(env, " must contain at least one value.", call. = FALSE)
  }
  if (isTRUE(numeric)) {
    values <- as.numeric(values)
    if (any(!is.finite(values))) {
      stop(env, " must contain finite numeric values.", call. = FALSE)
    }
  }

  values
}

# Prefer source loading during development, but allow installed-package use from
# a fresh clone or external validation directory.
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

# Extract the garbage collector's current maximum memory column in MB. This is a
# coarse signal, not a full profiler, but is useful across repeated runs.
gc_max_mb <- function() {
  gc_result <- gc()
  sum(gc_result[, 7], na.rm = TRUE)
}

# Summarize the returned POPMAPS list without converting every layer to files.
summarize_result <- function(result) {
  data.frame(
    output_layers = length(result),
    result_object_mb = as.numeric(utils::object.size(result)) / 1024^2,
    finite_output_values = sum(vapply(result, function(layer) {
      sum(is.finite(layer))
    }, numeric(1))),
    stringsAsFactors = FALSE
  )
}

# Resolve inputs, outputs, and fixed model settings before reading rasters.
args <- commandArgs(trailingOnly = TRUE)
raster_path <- read_arg_or_env(args, 1, "POPMAPS_BENCH_RASTER", "ASLO raster path")
locations_path <- read_arg_or_env(args, 2, "POPMAPS_BENCH_LOCS", "ASLO locations path")
output_dir <- read_arg_or_env(
  args,
  3,
  "POPMAPS_BENCH_OUTPUT_DIR",
  "output directory",
  required = FALSE,
  default = file.path(tempdir(), "popmaps2-aslo-benchmark")
)

load_popmaps2()
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("The terra package is required for ASLO benchmarking.", call. = FALSE)
}

raster_path <- normalizePath(raster_path, mustWork = TRUE)
locations_path <- normalizePath(locations_path, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

aggregate_facts <- as.integer(read_csv_values("POPMAPS_BENCH_AGGREGATES", "16,4,1", numeric = TRUE))
if (any(aggregate_facts < 1)) {
  stop("POPMAPS_BENCH_AGGREGATES must contain positive whole numbers.", call. = FALSE)
}
surfaces <- read_csv_values("POPMAPS_BENCH_SURFACES", "G,C")
if (any(!surfaces %in% c("G", "C"))) {
  stop("POPMAPS_BENCH_SURFACES must contain only G and/or C.", call. = FALSE)
}

empirical_pt_dist <- read_numeric_env("POPMAPS_BENCH_EMPIRICAL_PT_DIST", 5)
num_sites <- read_integer_env("POPMAPS_BENCH_NUM_SITES", 15)
num_tested <- read_integer_env("POPMAPS_BENCH_NUM_TESTED", 4)
popmod <- read_numeric_env("POPMAPS_BENCH_POPMOD", -0.05)
threshold <- read_numeric_env("POPMAPS_BENCH_THRESHOLD", 0)

# Read the base raster once and aggregate inside the loop. This makes the
# benchmark rows comparable and avoids repeated disk reads for each surface.
message("Reading ASLO raster: ", raster_path)
base_raster <- terra::rast(raster_path)
message("Reading ASLO locations: ", locations_path)
locations <- utils::read.table(
  locations_path,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

rows <- list()
row_idx <- 1L
for (aggregate_fact in aggregate_facts) {
  # Each aggregation factor represents a different problem size. The full
  # raster can be included with aggregate factor 1 when the machine can handle
  # it.
  raster_surface <- base_raster
  if (aggregate_fact > 1) {
    message("Aggregating raster by factor ", aggregate_fact)
    raster_surface <- terra::aggregate(raster_surface, fact = aggregate_fact, fun = mean, na.rm = TRUE)
  }

  raster_rows <- terra::nrow(raster_surface)
  raster_cols <- terra::ncol(raster_surface)
  raster_cells <- terra::ncell(raster_surface)
  raster_non_na <- sum(!is.na(terra::values(raster_surface, mat = FALSE)))

  for (surface in surfaces) {
    # Run each requested surface on the same raster/problem size. Errors are
    # captured into the output table so long benchmark sweeps can continue.
    message("Benchmarking aggregate=", aggregate_fact, ", surface=", surface)
    gc(reset = TRUE)
    result <- NULL
    error_message <- NA_character_
    runtime <- system.time({
      result <- tryCatch(
        popmaps2::popmaps(
          input_raster = raster_surface,
          input_locs = locations,
          surface = surface,
          empirical_pt_dist = empirical_pt_dist,
          num_sites = num_sites,
          num_tested = num_tested,
          popmod = popmod,
          threshold = threshold
        ),
        error = function(err) {
          error_message <<- conditionMessage(err)
          NULL
        }
      )
    })
    max_gc_mb <- gc_max_mb()

    result_summary <- if (is.null(result)) {
      data.frame(output_layers = NA_integer_, result_object_mb = NA_real_, finite_output_values = NA_real_)
    } else {
      summarize_result(result)
    }

    rows[[row_idx]] <- cbind(
      data.frame(
        created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        raster_path = raster_path,
        locations_path = locations_path,
        aggregate_fact = aggregate_fact,
        surface = surface,
        rows = raster_rows,
        cols = raster_cols,
        cells = raster_cells,
        non_na_cells = raster_non_na,
        empirical_pt_dist = empirical_pt_dist,
        num_sites = num_sites,
        num_tested = num_tested,
        popmod = popmod,
        threshold = threshold,
        elapsed_seconds = unname(runtime[["elapsed"]]),
        user_seconds = unname(runtime[["user.self"]]),
        system_seconds = unname(runtime[["sys.self"]]),
        gc_max_mb = max_gc_mb,
        error = error_message,
        stringsAsFactors = FALSE
      ),
      result_summary,
      popmaps_resource_row(resource_config)
    )
    row_idx <- row_idx + 1L
  }
}

benchmark <- do.call(rbind, rows)
output_path <- file.path(output_dir, paste0("aslo-benchmark-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"))
utils::write.csv(benchmark, output_path, row.names = FALSE)
print(benchmark)
message("Wrote benchmark summary: ", output_path)
