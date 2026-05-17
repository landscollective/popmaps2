#!/usr/bin/env Rscript

truthy_env <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

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

summarize_layers <- function(raster_output) {
  rows <- vector("list", terra::nlyr(raster_output))

  for (idx in seq_len(terra::nlyr(raster_output))) {
    values <- terra::values(raster_output[[idx]], mat = FALSE)
    finite_values <- values[is.finite(values)]

    rows[[idx]] <- data.frame(
      layer = names(raster_output)[idx],
      rows = terra::nrow(raster_output),
      cols = terra::ncol(raster_output),
      cells = terra::ncell(raster_output),
      non_na = sum(!is.na(values)),
      min = if (length(finite_values)) min(finite_values) else NA_real_,
      max = if (length(finite_values)) max(finite_values) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

args <- commandArgs(trailingOnly = TRUE)

raster_path <- read_arg_or_env(args, 1, "POPMAPS_ASLO_RASTER", "ASLO raster path")
locations_path <- read_arg_or_env(args, 2, "POPMAPS_ASLO_LOCS", "ASLO locations path")
output_dir <- read_arg_or_env(
  args,
  3,
  "POPMAPS_ASLO_OUTPUT_DIR",
  "output directory",
  required = FALSE,
  default = file.path(tempdir(), "popmaps2-aslo-validation")
)

if (!requireNamespace("popmaps2", quietly = TRUE)) {
  stop(
    "The popmaps2 package must be installed before running this script. ",
    "From the repository root, run: R CMD INSTALL .",
    call. = FALSE
  )
}
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("The terra package is required for ASLO validation.", call. = FALSE)
}

raster_path <- normalizePath(raster_path, mustWork = TRUE)
locations_path <- normalizePath(locations_path, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

aggregate_fact <- read_integer_env("POPMAPS_ASLO_AGGREGATE", 1)
surface <- Sys.getenv("POPMAPS_ASLO_SURFACE", unset = "G")
empirical_pt_dist <- read_numeric_env("POPMAPS_ASLO_EMPIRICAL_PT_DIST", 5)
num_sites <- read_integer_env("POPMAPS_ASLO_NUM_SITES", 15)
num_tested <- read_integer_env("POPMAPS_ASLO_NUM_TESTED", 4)
popmod <- read_numeric_env("POPMAPS_ASLO_POPMOD", -0.05)
threshold <- read_numeric_env("POPMAPS_ASLO_THRESHOLD", 0)
ncore <- read_integer_env("POPMAPS_ASLO_NCORE", 1)
write_rasters <- truthy_env(Sys.getenv("POPMAPS_ASLO_WRITE_RASTERS", unset = "false"))
save_rds <- truthy_env(Sys.getenv("POPMAPS_ASLO_SAVE_RDS", unset = "false"))

message("Reading ASLO raster: ", raster_path)
raster_surface <- terra::rast(raster_path)
if (aggregate_fact > 1) {
  message("Aggregating raster by factor ", aggregate_fact)
  raster_surface <- terra::aggregate(raster_surface, fact = aggregate_fact, fun = mean, na.rm = TRUE)
}

message("Reading ASLO locations: ", locations_path)
locations <- utils::read.table(
  locations_path,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

message("Running popmaps2 validation")
runtime <- system.time({
  result <- popmaps2::popmaps(
    input_raster = raster_surface,
    input_locs = locations,
    surface = surface,
    empirical_pt_dist = empirical_pt_dist,
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    threshold = threshold,
    ncore = ncore
  )
})

raster_output <- popmaps2::popmaps_rast(result, raster_surface)
layer_summary <- summarize_layers(raster_output)

run_summary <- data.frame(
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  raster_path = raster_path,
  locations_path = locations_path,
  aggregate_fact = aggregate_fact,
  surface = surface,
  empirical_pt_dist = empirical_pt_dist,
  num_sites = num_sites,
  num_tested = num_tested,
  popmod = popmod,
  threshold = threshold,
  ncore = ncore,
  elapsed_seconds = unname(runtime[["elapsed"]]),
  user_seconds = unname(runtime[["user.self"]]),
  system_seconds = unname(runtime[["sys.self"]]),
  stringsAsFactors = FALSE
)

run_summary_path <- file.path(output_dir, "aslo-run-summary.csv")
layer_summary_path <- file.path(output_dir, "aslo-layer-summary.csv")
utils::write.csv(run_summary, run_summary_path, row.names = FALSE)
utils::write.csv(layer_summary, layer_summary_path, row.names = FALSE)

manifest <- NULL
if (write_rasters) {
  raster_dir <- file.path(output_dir, "rasters")
  manifest <- popmaps2::write_popmaps(
    pop_raster_list = result,
    input_raster = raster_surface,
    dir = raster_dir,
    prefix = "aslo",
    overwrite = TRUE
  )
  utils::write.csv(manifest, file.path(output_dir, "aslo-raster-manifest.csv"), row.names = FALSE)
}

if (save_rds) {
  saveRDS(
    list(
      run_summary = run_summary,
      layer_summary = layer_summary,
      result = result
    ),
    file.path(output_dir, "aslo-validation-result.rds"),
    version = 2
  )
}

print(run_summary)
print(layer_summary)
message("Wrote run summary: ", run_summary_path)
message("Wrote layer summary: ", layer_summary_path)
if (!is.null(manifest)) {
  message("Wrote raster manifest: ", file.path(output_dir, "aslo-raster-manifest.csv"))
}
