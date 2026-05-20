#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE))
} else {
  file.path(getwd(), "tools")
}
source(file.path(script_dir, "popmaps-script-utils.R"))
resource_config <- popmaps_configure_script_resources("POPMAPS_EXAMPLE_SURFACE")

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

args <- commandArgs(trailingOnly = TRUE)
output_dir <- read_arg_or_env(
  args,
  1,
  "POPMAPS_EXAMPLE_SURFACE_OUTPUT_DIR",
  required = FALSE,
  default = file.path(getwd(), "local_validation", "example_surface_comparison")
)
aggregate_fact <- read_integer_env("POPMAPS_EXAMPLE_SURFACE_AGGREGATE", 120)
n_blocks <- read_integer_env("POPMAPS_EXAMPLE_SURFACE_N_BLOCKS", 3)
spatial_block_repeats <- read_integer_env("POPMAPS_EXAMPLE_SURFACE_BLOCK_REPEATS", 2)
spatial_block_seed <- read_integer_env("POPMAPS_EXAMPLE_SURFACE_BLOCK_SEED", 1, allow_zero = TRUE)

load_popmaps2()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

data("hija_raster", package = "popmaps2")
data("hija_struc", package = "popmaps2")

example_raster <- raster::aggregate(hija_raster, fact = aggregate_fact)
resolution <- raster::res(example_raster)[[1]]

geographic_surface <- popmaps2::surface_from_points(
  points = hija_struc,
  resolution = resolution,
  buffer = resolution,
  surface = "G"
)
suitability_surface <- popmaps2::prepare_popmaps_surface(
  input_raster = example_raster,
  surface = "C",
  surface_values = "suitability"
)

comparison <- popmaps2::compare_popmaps_surfaces(
  input_locs = hija_struc,
  surfaces = list(
    geographic_from_points = geographic_surface,
    sdm_suitability = suitability_surface
  ),
  validation = "spatial_block",
  n_blocks = n_blocks,
  spatial_block_repeats = spatial_block_repeats,
  spatial_block_seed = spatial_block_seed,
  num_sites = 5,
  num_tested = 2,
  near_best_tolerance = 0.05,
  quiet = TRUE
)

manifest <- popmaps2::write_surface_comparison_report(
  comparison,
  dir = output_dir,
  prefix = "hija-example-surface-comparison",
  overwrite = TRUE
)

run_summary <- cbind(
  data.frame(
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    output_dir = output_dir,
    aggregate_fact = aggregate_fact,
    n_blocks = n_blocks,
    spatial_block_repeats = spatial_block_repeats,
    spatial_block_seed = spatial_block_seed,
    report = manifest$report,
    stringsAsFactors = FALSE
  ),
  popmaps_resource_row(resource_config)
)
utils::write.csv(
  run_summary,
  file.path(output_dir, "hija-example-surface-comparison-run-summary.csv"),
  row.names = FALSE
)

print(comparison)
message("Wrote example report: ", manifest$report)
