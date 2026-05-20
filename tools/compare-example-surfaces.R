#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE))
} else {
  file.path(getwd(), "tools")
}
source(file.path(script_dir, "popmaps-script-utils.R"))
resource_config <- popmaps_configure_script_resources("POPMAPS_SURFACE")

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

read_numeric_env <- function(env, default, lower = -Inf) {
  value <- as.numeric(Sys.getenv(env, unset = as.character(default)))
  if (length(value) != 1 || !is.finite(value) || value < lower) {
    stop(env, " must be one finite number >= ", lower, ".", call. = FALSE)
  }

  value
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

read_species_env <- function(env) {
  species <- strsplit(Sys.getenv(env, unset = ""), ",", fixed = TRUE)[[1]]
  species <- trimws(species)
  species[nzchar(species)]
}

find_example_pairs <- function(input_dir, species_filter = character()) {
  rasters <- list.files(input_dir, pattern = "_avg[.]asc$", full.names = TRUE)
  rows <- lapply(rasters, function(raster_path) {
    species <- sub("_avg[.]asc$", "", basename(raster_path))
    if (length(species_filter) > 0 && !species %in% species_filter) {
      return(NULL)
    }
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

  pairs <- do.call(rbind, rows)
  pairs[order(pairs$species), , drop = FALSE]
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

markdown_table <- function(x, columns) {
  x <- x[, columns, drop = FALSE]
  x[] <- lapply(x, function(value) {
    if (is.numeric(value)) {
      signif(value, 4)
    } else {
      value
    }
  })

  header <- paste("|", paste(names(x), collapse = " | "), "|")
  divider <- paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|")
  rows <- apply(x, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  c(header, divider, rows)
}

default_num_sites <- function(n_locations, max_num_sites) {
  n_training <- n_locations - 1L
  candidates <- unique(c(5L, round(n_training * 0.5), min(max_num_sites, n_training), n_training))
  sort(candidates[candidates >= 2L & candidates <= n_training])
}

default_num_tested <- function(num_sites, max_num_tested) {
  seq.int(2L, min(max_num_tested, min(num_sites)))
}

flatten_grid <- function(grid, species, validation, surface_name) {
  data.frame(
    species = species,
    validation = validation,
    surface_name = surface_name,
    surface = grid$surface,
    surface_values = grid$surface_values,
    distance_units = grid$distance_units,
    reference_distance = grid$reference_distance,
    empirical_pt_dist = paste(signif(grid$empirical_pt_dist, 5), collapse = "; "),
    popmod = paste(signif(grid$popmod, 5), collapse = "; "),
    num_sites = paste(grid$num_sites, collapse = "; "),
    num_tested = paste(grid$num_tested, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

flatten_report_manifest <- function(manifest, species, validation) {
  data.frame(
    species = species,
    validation = validation,
    report = manifest$report,
    summary = manifest$tables[["summary"]],
    support = manifest$tables[["support"]],
    grids = manifest$tables[["grids"]],
    figures = paste(manifest$figures, collapse = "; "),
    tuning_results = paste(manifest$tuning_results, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

plot_surface_scores <- function(summary, path) {
  species <- sort(unique(summary$species))
  validations <- unique(summary$validation)
  surfaces <- unique(summary$surface_name)
  colors <- stats::setNames(grDevices::hcl.colors(length(surfaces), "Dark 3"), surfaces)

  png(path, width = 1500, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, length(validations)), mar = c(6, 5, 4, 1), oma = c(0, 0, 0, 0))
  for (validation in validations) {
    rows <- summary[summary$validation == validation, , drop = FALSE]
    score_range <- range(rows$score, finite = TRUE)
    graphics::plot(
      NA,
      xlim = c(0.5, length(species) + 0.5),
      ylim = score_range,
      xaxt = "n",
      xlab = "",
      ylab = paste("Validation", unique(rows$primary_metric)),
      main = validation
    )
    graphics::axis(1, at = seq_along(species), labels = species, las = 2)
    offsets <- seq(-0.16, 0.16, length.out = length(surfaces))
    for (idx in seq_along(surfaces)) {
      surface_name <- surfaces[[idx]]
      surface_rows <- rows[rows$surface_name == surface_name, , drop = FALSE]
      x <- match(surface_rows$species, species) + offsets[[idx]]
      graphics::points(x, surface_rows$score, pch = 19, col = colors[[surface_name]], cex = 1.3)
      graphics::segments(x, score_range[1], x, surface_rows$score, col = grDevices::adjustcolor(colors[[surface_name]], 0.35))
    }
    if (validation == validations[[1]]) {
      graphics::legend("topright", legend = surfaces, pch = 19, col = colors, bty = "n", cex = 0.85)
    }
  }
}

plot_percent_from_best <- function(summary, path) {
  species <- sort(unique(summary$species))
  validations <- unique(summary$validation)
  surfaces <- unique(summary$surface_name)
  colors <- stats::setNames(grDevices::hcl.colors(length(surfaces), "Dark 3"), surfaces)

  png(path, width = 1500, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, length(validations)), mar = c(6, 5, 4, 1), oma = c(0, 0, 0, 0))
  for (validation in validations) {
    rows <- summary[summary$validation == validation, , drop = FALSE]
    y_max <- max(5, rows$percent_from_best, na.rm = TRUE)
    graphics::plot(
      NA,
      xlim = c(0.5, length(species) + 0.5),
      ylim = c(0, y_max),
      xaxt = "n",
      xlab = "",
      ylab = "Percent worse than best",
      main = validation
    )
    graphics::axis(1, at = seq_along(species), labels = species, las = 2)
    graphics::abline(h = 5, lty = 2, col = "gray45")
    offsets <- seq(-0.16, 0.16, length.out = length(surfaces))
    for (idx in seq_along(surfaces)) {
      surface_name <- surfaces[[idx]]
      surface_rows <- rows[rows$surface_name == surface_name, , drop = FALSE]
      x <- match(surface_rows$species, species) + offsets[[idx]]
      graphics::points(x, surface_rows$percent_from_best, pch = 19, col = colors[[surface_name]], cex = 1.3)
      graphics::segments(x, 0, x, surface_rows$percent_from_best, col = grDevices::adjustcolor(colors[[surface_name]], 0.35))
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)
default_input_dir <- file.path(dirname(getwd()), "popmaps_test_data")
input_dir <- read_arg_or_env(args, 1, "POPMAPS_SURFACE_EXAMPLE_DIR", default = default_input_dir)
output_dir <- read_arg_or_env(
  args,
  2,
  "POPMAPS_SURFACE_OUTPUT_DIR",
  required = FALSE,
  default = file.path(input_dir, "surface_comparison_outputs")
)

load_popmaps2()
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("The terra package is required for empirical surface comparison.", call. = FALSE)
}

input_dir <- normalizePath(input_dir, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

species_filter <- read_species_env("POPMAPS_SURFACE_SPECIES")
validation_modes <- read_modes_env("POPMAPS_SURFACE_VALIDATION", "spatial_block")
aggregate_fact <- read_integer_env("POPMAPS_SURFACE_AGGREGATE", 8)
n_blocks <- read_integer_env("POPMAPS_SURFACE_N_BLOCKS", 4)
spatial_block_repeats <- read_integer_env("POPMAPS_SURFACE_BLOCK_REPEATS", 2)
spatial_block_seed <- read_integer_env("POPMAPS_SURFACE_BLOCK_SEED", 1, allow_zero = TRUE)
near_best_tolerance <- read_numeric_env("POPMAPS_SURFACE_NEAR_BEST_TOLERANCE", 0.05, lower = 0)
max_num_sites <- read_integer_env("POPMAPS_SURFACE_MAX_NUM_SITES", 15)
max_num_tested <- read_integer_env("POPMAPS_SURFACE_MAX_NUM_TESTED", 5)
include_inverse <- truthy_env(Sys.getenv("POPMAPS_SURFACE_INCLUDE_INVERSE", unset = "false"))
write_surface_results <- truthy_env(Sys.getenv("POPMAPS_SURFACE_WRITE_RESULTS", unset = "true"))

pairs <- find_example_pairs(input_dir, species_filter = species_filter)
message("Found ", nrow(pairs), " empirical example datasets.")

run_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
summary_rows <- list()
support_rows <- list()
grid_rows <- list()
report_rows <- list()

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
  num_sites <- default_num_sites(nrow(locations), max_num_sites = max_num_sites)
  num_tested <- default_num_tested(num_sites, max_num_tested = max_num_tested)

  surfaces <- list(
    geographic = popmaps2::prepare_popmaps_surface(raster_surface, surface = "G"),
    sdm_suitability = popmaps2::prepare_popmaps_surface(
      raster_surface,
      surface = "C",
      surface_values = "suitability"
    )
  )
  if (include_inverse) {
    surfaces$sdm_as_resistance <- popmaps2::prepare_popmaps_surface(
      raster_surface,
      surface = "C",
      surface_values = "resistance"
    )
  }

  for (validation in validation_modes) {
    message("Comparing surfaces for ", species, " with validation = ", validation)
    comparison <- tryCatch(
      popmaps2::compare_popmaps_surfaces(
        input_locs = locations,
        surfaces = surfaces,
        num_sites = num_sites,
        num_tested = num_tested,
        validation = validation,
        n_blocks = n_blocks,
        spatial_block_repeats = if (validation == "spatial_block") spatial_block_repeats else 1,
        spatial_block_seed = if (validation == "spatial_block") spatial_block_seed else NULL,
        near_best_tolerance = near_best_tolerance,
        quiet = TRUE
      ),
      error = function(err) err
    )

    prefix <- paste(species, validation, "surface-comparison", run_stamp, sep = "-")
    if (inherits(comparison, "error")) {
      support <- data.frame(
        species = species,
        validation = validation,
        decision = "failed",
        best_surface = NA_character_,
        primary_metric = "rmse",
        best_score = NA_real_,
        near_best_tolerance = near_best_tolerance,
        n_near_best = NA_integer_,
        near_best_surfaces = NA_character_,
        error = conditionMessage(comparison),
        stringsAsFactors = FALSE
      )
      support_rows[[length(support_rows) + 1]] <- support
      report_rows[[length(report_rows) + 1]] <- data.frame(
        species = species,
        validation = validation,
        report = NA_character_,
        summary = NA_character_,
        support = file.path(output_dir, paste0(prefix, "-support.csv")),
        grids = NA_character_,
        figures = NA_character_,
        tuning_results = NA_character_,
        stringsAsFactors = FALSE
      )
      write_table(support, file.path(output_dir, paste0(prefix, "-support.csv")))
      next
    }

    comparison_report_dir <- file.path(output_dir, "comparisons", prefix)
    report_manifest <- popmaps2::write_surface_comparison_report(
      comparison = comparison,
      dir = comparison_report_dir,
      prefix = prefix,
      include_tuning_results = write_surface_results,
      overwrite = TRUE
    )
    report_rows[[length(report_rows) + 1]] <- flatten_report_manifest(
      report_manifest,
      species = species,
      validation = validation
    )

    summary <- cbind(
      data.frame(species = species, validation = validation, stringsAsFactors = FALSE),
      comparison$summary
    )
    support <- cbind(
      data.frame(species = species, validation = validation, stringsAsFactors = FALSE),
      comparison$support
    )
    grids <- do.call(rbind, lapply(names(comparison$grids), function(surface_name) {
      flatten_grid(comparison$grids[[surface_name]], species, validation, surface_name)
    }))

    summary_rows[[length(summary_rows) + 1]] <- summary
    support_rows[[length(support_rows) + 1]] <- support
    grid_rows[[length(grid_rows) + 1]] <- grids
  }
}

summary_table <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()
support_table <- if (length(support_rows) > 0) do.call(rbind, support_rows) else data.frame()
grid_table <- if (length(grid_rows) > 0) do.call(rbind, grid_rows) else data.frame()
report_table <- if (length(report_rows) > 0) do.call(rbind, report_rows) else data.frame()

summary_path <- file.path(output_dir, paste0("empirical-surface-comparison-summary-", run_stamp, ".csv"))
support_path <- file.path(output_dir, paste0("empirical-surface-comparison-support-", run_stamp, ".csv"))
grid_path <- file.path(output_dir, paste0("empirical-surface-comparison-grids-", run_stamp, ".csv"))
report_manifest_path <- file.path(output_dir, paste0("empirical-surface-comparison-reports-", run_stamp, ".csv"))
resource_path <- file.path(output_dir, paste0("empirical-surface-comparison-run-summary-", run_stamp, ".csv"))
write_table(summary_table, summary_path)
write_table(support_table, support_path)
write_table(grid_table, grid_path)
write_table(report_table, report_manifest_path)
write_table(
  cbind(
    data.frame(
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      input_dir = input_dir,
      output_dir = output_dir,
      aggregate_fact = aggregate_fact,
      validation_modes = paste(validation_modes, collapse = ","),
      n_datasets = nrow(pairs),
      include_inverse = include_inverse,
      stringsAsFactors = FALSE
    ),
    popmaps_resource_row(resource_config)
  ),
  resource_path
)

report_dir <- file.path(output_dir, paste0("report-", run_stamp))
figure_dir <- file.path(report_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(summary_table) > 0) {
  plot_surface_scores(summary_table, file.path(figure_dir, "surface-validation-scores.png"))
  plot_percent_from_best(summary_table, file.path(figure_dir, "surface-percent-from-best.png"))
}

report_path <- file.path(report_dir, "empirical-surface-comparison-report.md")
report <- c(
  "# Empirical Surface Comparison Report",
  "",
  paste0("Source summary table: `", summary_path, "`"),
  paste0("Source support table: `", support_path, "`"),
  paste0("Source grid table: `", grid_path, "`"),
  paste0("Source report manifest: `", report_manifest_path, "`"),
  paste0("Source run summary table: `", resource_path, "`"),
  "",
  paste0("Raster aggregation factor: ", aggregate_fact),
  paste0("Validation modes: ", paste(validation_modes, collapse = ", ")),
  paste0("Spatial block repeats: ", spatial_block_repeats),
  paste0(
    "Resource configuration: ",
    resource_config$threads,
    " of ",
    resource_config$available_threads,
    " detected logical processors; terra memfrac = ",
    resource_config$terra_memfrac
  ),
  "",
  "Lower RMSE values are better. `percent_from_best` is the percent increase in RMSE relative to the best-ranked surface for that species and validation design.",
  "`surfaces_indistinguishable` means at least two surfaces were within the near-best tolerance, so the data do not clearly support one surface over the other.",
  "Each completed species/validation comparison also writes the standard `write_surface_comparison_report()` CSVs, figures, and Markdown report in `comparisons/`.",
  "",
  "## Surface Support",
  "",
  if (nrow(support_table) > 0) {
    markdown_table(
      support_table,
      intersect(
        c("species", "validation", "decision", "best_surface", "best_score",
          "n_near_best", "near_best_surfaces", "error"),
        names(support_table)
      )
    )
  } else {
    "No completed comparisons."
  },
  "",
  "## Surface Ranking",
  "",
  if (nrow(summary_table) > 0) {
    markdown_table(
      summary_table,
      intersect(
        c("species", "validation", "rank", "surface_name", "surface",
          "surface_values", "score", "percent_from_best", "distance_units",
          "num_sites", "num_tested", "popmod", "empirical_pt_dist"),
        names(summary_table)
      )
    )
  } else {
    "No completed surface rankings."
  },
  "",
  "## Figures",
  "",
  "- `figures/surface-validation-scores.png`: validation RMSE by species and surface.",
  "- `figures/surface-percent-from-best.png`: relative support gap between candidate surfaces.",
  "",
  "## Per-Comparison Reports",
  "",
  if (nrow(report_table) > 0) {
    markdown_table(
      report_table,
      intersect(c("species", "validation", "report"), names(report_table))
    )
  } else {
    "No comparison reports were written."
  },
  ""
)
writeLines(report, report_path)
message("Wrote ", report_path)

if (nrow(support_table) > 0) {
  print(support_table)
}
