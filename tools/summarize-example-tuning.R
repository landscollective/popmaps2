#!/usr/bin/env Rscript

# Summarize empirical tuning outputs.
#
# Purpose:
# - read the latest CSV bundle produced by `validate-example-tuning.R`;
# - create cross-species summary tables;
# - draw diagnostic figures for validation score, near-best support, decay
#   scales, and parameter effects;
# - write a compact Markdown report for interpretation.
#
# Usage:
#   Rscript tools/summarize-example-tuning.R tuning_output_dir report_dir
#
# This script does not rerun models. It only summarizes existing tuning CSVs.

# Resolve shared utilities from this script location when available.
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE))
} else {
  file.path(getwd(), "tools")
}
source(file.path(script_dir, "popmaps-script-utils.R"))
resource_config <- popmaps_configure_script_resources("POPMAPS_TUNING")

# Positional argument, environment variable, default: in that order.
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

# Pick the newest matching CSV so users can point at a tuning output directory
# without manually copying timestamped filenames.
latest_file <- function(output_dir, pattern) {
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(files) < 1) {
    stop("No files matching '", pattern, "' found in ", output_dir, ".", call. = FALSE)
  }

  files[which.max(file.info(files)$mtime)]
}

# Use stable CSV options across all report inputs.
read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# CSV writer that also logs artifact paths during report generation.
write_table <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  message("Wrote ", path)
}

# Translate the breadth of near-best support into a compact interpretation.
support_label <- function(n_near_best, near_best_fraction) {
  ifelse(
    n_near_best <= 1 | near_best_fraction <= 0.10,
    "sharp",
    ifelse(near_best_fraction <= 0.25, "moderate", "broad")
  )
}

# Combine performance improvement and support breadth into a rough signal label.
tuning_signal <- function(best_vs_median_percent, support) {
  ifelse(
    best_vs_median_percent >= 25 & support == "sharp",
    "strong",
    ifelse(best_vs_median_percent >= 10 | support != "broad", "moderate", "weak")
  )
}

# Small Markdown table renderer for the final report.
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

# Keep species ordering consistent across plots.
ordered_species <- function(x) {
  sort(unique(x$species))
}

# Fixed validation colors make reports easier to compare across repeated runs.
validation_colors <- function(validation) {
  colors <- c(loo = "#1f77b4", spatial_block = "#d95f02")
  stats::setNames(colors[validation], validation)
}

# Plot best validation score by species. When repeated spatial blocks are
# available, show repeat-level standard deviation as error bars.
plot_best_score <- function(summary, path) {
  species <- ordered_species(summary)
  validations <- unique(summary$validation)
  colors <- validation_colors(validations)
  offsets <- seq(-0.16, 0.16, length.out = length(validations))
  metric <- unique(summary$primary_metric)
  y_label <- if (length(metric) == 1) paste("Best", metric) else "Best score"
  score_interval <- c(
    summary$best_score - summary$best_score_repeat_sd,
    summary$best_score + summary$best_score_repeat_sd
  )

  png(path, width = 1400, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(6, 5, 4, 2))
  graphics::plot(
    NA,
    xlim = c(0.5, length(species) + 0.5),
    ylim = range(c(summary$best_score, summary$median_score, score_interval), finite = TRUE),
    xaxt = "n",
    xlab = "",
    ylab = y_label,
    main = "Best validation error by species"
  )
  graphics::axis(1, at = seq_along(species), labels = species, las = 2)
  for (idx in seq_along(validations)) {
    validation <- validations[idx]
    rows <- summary[summary$validation == validation, , drop = FALSE]
    x <- match(rows$species, species) + offsets[idx]
    graphics::points(x, rows$best_score, pch = 19, col = colors[[validation]], cex = 1.4)
    graphics::segments(x, rows$median_score, x, rows$best_score, col = grDevices::adjustcolor(colors[[validation]], alpha.f = 0.45))
    has_sd <- is.finite(rows$best_score_repeat_sd)
    if (any(has_sd)) {
      graphics::arrows(
        x[has_sd],
        rows$best_score[has_sd] - rows$best_score_repeat_sd[has_sd],
        x[has_sd],
        rows$best_score[has_sd] + rows$best_score_repeat_sd[has_sd],
        angle = 90,
        code = 3,
        length = 0.04,
        col = colors[[validation]]
      )
      graphics::points(x[has_sd], rows$best_score[has_sd], pch = 19, col = colors[[validation]], cex = 1.4)
    }
  }
  graphics::legend("topright", legend = validations, pch = 19, col = colors, bty = "n")
}

# Plot the fraction of evaluated parameter combinations within the near-best
# tolerance; smaller fractions indicate sharper parameter support.
plot_near_best <- function(summary, path) {
  species <- ordered_species(summary)
  validations <- unique(summary$validation)
  colors <- validation_colors(validations)
  values <- vapply(validations, function(validation) {
    rows <- summary[summary$validation == validation, , drop = FALSE]
    rows$near_best_fraction[match(species, rows$species)]
  }, numeric(length(species)))
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1)
    colnames(values) <- validations
  }
  rownames(values) <- species

  png(path, width = 1400, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(6, 5, 4, 2))
  graphics::barplot(
    t(values),
    beside = TRUE,
    col = colors,
    border = NA,
    ylim = c(0, max(1, values, na.rm = TRUE)),
    ylab = "Near-best combinations / evaluated combinations",
    main = "Breadth of near-best support",
    las = 2,
    names.arg = species
  )
  graphics::abline(h = c(0.10, 0.25), lty = c(2, 3), col = "gray45")
  graphics::legend("topright", legend = validations, fill = colors, bty = "n")
}

# Plot interpretable distance-decay scales for the best parameter combination.
plot_decay <- function(summary, path) {
  decay <- summary[is.finite(summary$half_distance_km) & is.finite(summary$ten_pct_distance_km), , drop = FALSE]
  species <- ordered_species(decay)
  validations <- unique(decay$validation)
  colors <- validation_colors(validations)
  offsets <- seq(-0.16, 0.16, length.out = length(validations))

  png(path, width = 1400, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(6, 5, 4, 2))
  graphics::plot(
    NA,
    xlim = c(0.5, length(species) + 0.5),
    ylim = range(c(decay$half_distance_km, decay$ten_pct_distance_km), finite = TRUE),
    log = "y",
    xaxt = "n",
    xlab = "",
    ylab = "Distance decay scale (km, log scale)",
    main = "Best distance-decay scales"
  )
  graphics::axis(1, at = seq_along(species), labels = species, las = 2)
  for (idx in seq_along(validations)) {
    validation <- validations[idx]
    rows <- decay[decay$validation == validation, , drop = FALSE]
    x <- match(rows$species, species) + offsets[idx]
    graphics::segments(x, rows$half_distance_km, x, rows$ten_pct_distance_km, col = colors[[validation]], lwd = 2)
    graphics::points(x, rows$half_distance_km, pch = 19, col = colors[[validation]], cex = 1.2)
    graphics::points(x, rows$ten_pct_distance_km, pch = 17, col = colors[[validation]], cex = 1.2)
  }
  graphics::legend(
    "topleft",
    legend = c(paste(validations, "validation"), "50% weight", "10% weight"),
    col = c(colors, "black", "black"),
    pch = c(rep(19, length(validations)), 19, 17),
    bty = "n"
  )
}

# Plot mean loss from the best score across each parameter value. This helps
# diagnose whether tuning is sensitive to parameter changes.
plot_parameter_effects <- function(effects, path) {
  effects <- effects[is.finite(effects$value) & is.finite(effects$loss_from_best), , drop = FALSE]
  parameters <- unique(effects$parameter)
  species <- ordered_species(effects)
  species_colors <- grDevices::hcl.colors(length(species), "Dark 3")
  names(species_colors) <- species
  validation_lty <- c(loo = 1, spatial_block = 2)

  png(path, width = 1500, height = 1100, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(5, 5, 4, 1), oma = c(0, 0, 1, 0))
  for (parameter in parameters) {
    rows <- effects[effects$parameter == parameter, , drop = FALSE]
    graphics::plot(
      NA,
      xlim = range(rows$value, finite = TRUE),
      ylim = range(rows$loss_from_best, finite = TRUE),
      xlab = parameter,
      ylab = "Mean score - best score",
      main = parameter
    )
    if (parameter == parameters[[1]]) {
      graphics::legend(
        "topleft",
        legend = names(validation_lty),
        col = "gray30",
        lty = validation_lty,
        bty = "n",
        cex = 0.85
      )
    }
    groups <- unique(paste(rows$species, rows$validation, sep = "|"))
    for (group in groups) {
      parts <- strsplit(group, "|", fixed = TRUE)[[1]]
      species_name <- parts[[1]]
      validation <- parts[[2]]
      group_rows <- rows[rows$species == species_name & rows$validation == validation, , drop = FALSE]
      group_rows <- group_rows[order(group_rows$value), , drop = FALSE]
      graphics::lines(
        group_rows$value,
        group_rows$loss_from_best,
        type = "b",
        pch = 19,
        col = grDevices::adjustcolor(species_colors[[species_name]], alpha.f = 0.65),
        lty = validation_lty[[validation]]
      )
    }
  }
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  graphics::legend("bottom", legend = species, col = species_colors, pch = 19, horiz = TRUE, bty = "n", cex = 0.8)
}

# Resolve the output bundle to summarize and create the report directory.
args <- commandArgs(trailingOnly = TRUE)
output_dir <- read_arg_or_env(args, 1, "POPMAPS_TUNING_OUTPUT_DIR")
report_dir <- read_arg_or_env(
  args,
  2,
  "POPMAPS_TUNING_REPORT_DIR",
  required = FALSE,
  default = file.path(output_dir, paste0("report-", format(Sys.time(), "%Y%m%d-%H%M%S")))
)

output_dir <- normalizePath(output_dir, mustWork = TRUE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report_dir <- normalizePath(report_dir, mustWork = TRUE)
figure_dir <- file.path(report_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

best_path <- latest_file(output_dir, "^empirical-tuning-best-.*[.]csv$")
overview_path <- latest_file(output_dir, "^empirical-tuning-overview-.*[.]csv$")
range_path <- latest_file(output_dir, "^empirical-tuning-parameter-ranges-.*[.]csv$")
effect_path <- latest_file(output_dir, "^empirical-tuning-parameter-effects-.*[.]csv$")

best <- read_csv(best_path)
overview <- read_csv(overview_path)
parameter_ranges <- read_csv(range_path)
parameter_effects <- read_csv(effect_path)

summary <- merge(
  best,
  overview[, c("species", "search", "validation", "primary_metric", "metric_goal",
               "n_combinations", "n_evaluated", "n_complete", "best_score",
               "median_score", "worst_score", "best_vs_median_percent",
               "best_vs_worst_percent", "near_best_tolerance", "n_near_best")],
  by = c("species", "search", "validation"),
  suffixes = c("", "_overview")
)
# Add derived interpretation fields used in both CSV and Markdown outputs.
summary$near_best_fraction <- summary$n_near_best / summary$n_evaluated
summary$best_score_repeat_sd <- vapply(seq_len(nrow(summary)), function(row_idx) {
  metric_sd_col <- paste0(summary$primary_metric[row_idx], "_repeat_sd")
  if (metric_sd_col %in% names(summary)) {
    summary[[metric_sd_col]][row_idx]
  } else {
    NA_real_
  }
}, numeric(1))
summary$support <- support_label(summary$n_near_best, summary$near_best_fraction)
summary$tuning_signal <- tuning_signal(summary$best_vs_median_percent, summary$support)
summary <- summary[order(summary$validation, summary$species), , drop = FALSE]

diversity_rows <- lapply(split(summary, summary$validation), function(rows) {
  # Count how many distinct best-parameter signatures appear across species.
  # If every species chooses the same parameters, tuning may be less informative
  # than expected; diversity here supports species-specific tuning.
  signatures <- paste(
    rows$num_sites,
    rows$num_tested,
    signif(rows$popmod, 6),
    signif(rows$empirical_pt_dist, 6),
    sep = "|"
  )
  data.frame(
    validation = rows$validation[[1]],
    n_species = length(unique(rows$species)),
    unique_best_parameter_sets = length(unique(signatures)),
    unique_num_sites = length(unique(rows$num_sites)),
    unique_num_tested = length(unique(rows$num_tested)),
    unique_popmod = length(unique(signif(rows$popmod, 6))),
    unique_empirical_pt_dist = length(unique(signif(rows$empirical_pt_dist, 6))),
    stringsAsFactors = FALSE
  )
})
parameter_diversity <- do.call(rbind, diversity_rows)
rownames(parameter_diversity) <- NULL

# Write the cleaned summary tables before plotting so partial report generation
# still leaves useful CSV outputs if plotting fails.
write_table(summary, file.path(report_dir, "empirical-tuning-summary.csv"))
write_table(parameter_diversity, file.path(report_dir, "empirical-tuning-parameter-diversity.csv"))
write_table(parameter_ranges, file.path(report_dir, "empirical-tuning-parameter-ranges.csv"))
write_table(
  cbind(
    data.frame(
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      output_dir = output_dir,
      report_dir = report_dir,
      stringsAsFactors = FALSE
    ),
    popmaps_resource_row(resource_config)
  ),
  file.path(report_dir, "empirical-tuning-report-run-summary.csv")
)
write_table(parameter_effects, file.path(report_dir, "empirical-tuning-parameter-effects.csv"))

plot_best_score(summary, file.path(figure_dir, "best-validation-score.png"))
plot_near_best(summary, file.path(figure_dir, "near-best-support.png"))
plot_decay(summary, file.path(figure_dir, "distance-decay-scales.png"))
plot_parameter_effects(parameter_effects, file.path(figure_dir, "parameter-effects.png"))

summary_columns <- intersect(
  c("species", "validation", "n_validation_repeats", "best_score",
    "best_score_repeat_sd", "near_best_fraction", "support",
    "tuning_signal", "num_sites", "num_tested", "popmod",
    "half_distance_km", "ten_pct_distance_km", "empirical_pt_dist"),
  names(summary)
)
report_path <- file.path(report_dir, "empirical-tuning-report.md")
report <- c(
  "# Empirical Tuning Report",
  "",
  paste0("Source best table: `", best_path, "`"),
  paste0("Source overview table: `", overview_path, "`"),
  paste0("Source parameter ranges table: `", range_path, "`"),
  paste0("Source parameter effects table: `", effect_path, "`"),
  "",
  "Lower RMSE/MAE/Hellinger values are better; dominant accuracy and dominant probability are maximized when used as the primary metric.",
  "`near_best_fraction` is the fraction of evaluated parameter combinations within the near-best tolerance. Smaller values mean sharper tuning support.",
  "",
  "## Species Summary",
  "",
  markdown_table(summary, summary_columns),
  "",
  "## Parameter Diversity",
  "",
  markdown_table(parameter_diversity, names(parameter_diversity)),
  "",
  "## Figures",
  "",
  "- `figures/best-validation-score.png`: best validation score by species and validation design. Error bars show repeat-level standard deviation when repeated spatial blocks are available.",
  "- `figures/near-best-support.png`: whether each species has sharp or broad near-best support.",
  "- `figures/distance-decay-scales.png`: best 50% and 10% distance-decay scales.",
  "- `figures/parameter-effects.png`: average loss from the best score across each tuning parameter.",
  ""
)
writeLines(report, report_path)
message("Wrote ", report_path)
print_columns <- intersect(
  c("species", "validation", "n_validation_repeats", "best_score",
    "best_score_repeat_sd", "near_best_fraction", "support", "tuning_signal"),
  names(summary)
)
print(summary[, print_columns, drop = FALSE])
