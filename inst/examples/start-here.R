# This script is the smallest complete popmaps2 workflow:
# 1. prepare empirical ancestry locations and a prediction raster,
# 2. tune a small parameter grid,
# 3. run POPMAPS with the selected parameters,
# 4. write rasters and a map figure.
#
# It uses the small bundled Hesperidanthus jaegeri example data so it can run
# quickly after installation.

library(popmaps2)

# The bundled raster is intentionally aggregated here so the example finishes
# quickly. For real analyses, use the raster resolution appropriate for the
# management question and available computing resources.
example_raster <- raster::aggregate(hija_raster, fact = 240)
example_locs <- hija_struc

# Tune a tiny, transparent grid first. Larger analyses should use
# suggest_tuning_grid(), adaptive_tune_popmaps(), or a biologically motivated
# grid that covers plausible dispersal/gene-flow scales.
tuning <- tune_popmaps(
  input_raster = example_raster,
  input_locs = example_locs,
  surface = "G",
  empirical_pt_dist = c(0, 5),
  num_sites = c(5, 8),
  num_tested = c(2, 3),
  popmod = c(-0.01, -0.05),
  threshold = 0,
  quiet = FALSE
)

print(tuning)

# Run the final map with the best parameter combination selected by RMSE.
# Lower RMSE is better because it measures prediction error in withheld
# empirical ancestry coefficients.
best <- tuning$best
aps <- popmaps(
  input_raster = example_raster,
  input_locs = example_locs,
  surface = "G",
  empirical_pt_dist = best$empirical_pt_dist,
  num_sites = best$num_sites,
  num_tested = best$num_tested,
  popmod = best$popmod,
  threshold = 0,
  quiet = FALSE
)

# Convert the legacy POPMAPS list output to a named terra SpatRaster for modern
# spatial workflows.
aps_raster <- popmaps_rast(aps, example_raster)
print(aps_raster)

# Write reproducible outputs to a temporary directory by default. Change this
# path to a project folder when running your own analysis.
output_dir <- file.path(tempdir(), "popmaps2-start-here")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

written_layers <- write_popmaps(
  pop_raster_list = aps,
  input_raster = example_raster,
  dir = output_dir,
  prefix = "hija",
  overwrite = TRUE
)
print(written_layers)

map_path <- file.path(output_dir, "hija-ancestry-map.png")
write_popmaps_plot(
  pop_raster_list = aps,
  input_raster = example_raster,
  path = map_path,
  input_locs = example_locs,
  style = "manuscript",
  type = "confidence",
  overwrite = TRUE
)

message("Example outputs written to: ", normalizePath(output_dir, mustWork = FALSE))
