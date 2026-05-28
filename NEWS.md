# popmaps2 0.1.0

- Prepared the first public alpha release of `popmaps2`. The package is
  installable, tested, documented, and suitable for validation work, but the API
  may still change before CRAN submission and the software manuscript.
- Added a manual full-platform R-CMD-check workflow for release-readiness checks
  across Ubuntu, macOS, Windows, and R-devel.
- Added a bundled start-here example script, early coordinate/raster/ancestry
  validation, clearer input errors, and optional progress messages for
  `popmaps()` and `tune_popmaps()`.

- Seeded package from the USGS POPMAPS 1.03 release.
- Renamed the package to `popmaps2`.
- Added package metadata for GitHub-based development.
- Split core dependencies from optional legacy modeling/plotting dependencies.
- Added clearer errors when optional packages are required.
- Fixed default argument matching in legacy plotting functions.
- Added README and roadmap documentation for modernization and public release.
- Added terra-first input preparation and validation for raster and location inputs.
- Allowed descriptive location column names while preserving the legacy column-order contract.
- Updated `anc_extract()` to use terra cell lookup and clearer coordinate validation.
- Removed the retired `rgeos` package from metadata and updated `popmap_viz()` boundary drawing.
- Added `validate_popmaps_baseline()` for repeatable comparison against a frozen POPMAPS 1.03 reference result.
- Added an optimized geographic-distance engine for `popmaps(surface = "G")` that precomputes cell-site and empirical-site distances while preserving POPMAPS 1.03 outputs on validation cases.
- Added `popmaps_rast()` and `write_popmaps()` to convert `popmaps()` output to `terra` rasters and export GeoTIFF layers.
- Added `tools/validate-aslo-local.R` for local, larger-scale ASLO validation without committing raw project data.
- Added `tune_popmaps()` for fast leave-one-site-out tuning of geographic-distance POPMAPS parameters with fold-level diagnostics and summary metrics.
- Added vectorized `empirical_pt_dist` tuning plus `suggest_tuning_grid()` and `adaptive_tune_popmaps()` for data-adaptive parameter exploration.
- Added `half_distance_km` and `ten_pct_distance_km` tuning summaries plus spatial-block validation for more biologically interpretable parameter selection.
- Added `diagnose_tuning()` to summarize tuning strength, near-best parameter support, and parameter effects.
- Added `tools/validate-example-tuning.R` for repeatable local tuning validation across empirical example datasets kept outside the package repository.
- Added `tools/summarize-example-tuning.R` to create empirical tuning summary tables, a markdown report, and diagnostic plots.
- Added repeated spatial-block validation with reproducible rotated spatial partitions and repeat-level uncertainty summaries.
- Documented the scientific contract for comparing geographic-distance (`surface = "G"`) and SDM suitability-weighted least-cost (`surface = "C"`) surfaces.
- Added `prepare_popmaps_surface()` to declare candidate raster semantics, including suitability, conductance, and resistance-to-conductance transformations.
- Added internal least-cost distance helpers that mirror POPMAPS 1.03 `gdistance` conductance distances on small validation rasters.
- Extended `tune_popmaps(surface = "C")` to score suitability-, conductance-, and resistance-weighted least-cost surfaces with the modern distance helper.
- Added `compare_popmaps_surfaces()` for matched predictive comparison of user-supplied candidate surfaces.
- Added `suggest_surface_tuning_grid()` and surface-specific defaults in `compare_popmaps_surfaces()` so geographic and least-cost surfaces can be tuned on their own distance scales.
- Added `plot_surface_comparison()` and `write_surface_comparison_report()` to turn candidate-surface comparisons into diagnostic figures, CSVs, and Markdown reports.
- Added `tools/compare-example-surfaces.R` for local empirical G-vs-SDM surface comparisons, summary tables, reports, and plots across uncommitted example datasets.
- Updated `tools/compare-example-surfaces.R` so each completed species/validation comparison writes the standard surface-comparison report artifacts.
- Replaced the full `popmaps(surface = "C")` least-cost interpolation engine with the internal distance helper, removing the `gdistance` requirement from the main ancestry-surface workflow while preserving the original suitability-as-conductance default.
- Added shared resource configuration for local validation/reporting scripts, including cross-platform processor detection, conservative thread defaults for system libraries, terra memory settings, and run-summary reporting.
- Added `tools/benchmark-aslo-local.R` for repeatable small-, medium-, and full-scale ASLO runtime and memory checks.
- Added input converter helpers: `surface_from_points()`, `locs_from_sf()`, `surfaces_from_raster_stack()`, `surface_from_eems()`, and `surface_from_feems()`.
- Added a surface-comparison vignette and a bundled-data example script that uses `surface_from_points()` plus `write_surface_comparison_report()`.
- Added a pkgdown scaffold, clearer vignette installation notes, and runnable examples for surface input converter help pages.
- Added a manual pkgdown deployment workflow, reduced routine R-CMD-check triggers to pull requests/manual runs, cleaned the README, and annotated local validation scripts.
- Corrected default `popmaps()` raster-cell coordinates to use true raster cell centers, with `legacy_compat = TRUE` retained for POPMAPS 1.03 reproducibility.
- Updated default `surface = "C"` interpolation to select candidate empirical sites by least-cost distance, while legacy geographic prefiltering remains available only through `legacy_compat = TRUE`.
- Added unit-neutral `half_distance` and `ten_pct_distance` tuning summaries while retaining the older `_km` columns for compatibility.
- Changed optional conductance rescaling to divide by maximum conductance so zero-valued barriers remain zero instead of turning the minimum positive conductance into a barrier.
- Added `plot_popmaps()` and `write_popmaps_plot()` as the modern POPMAPS visualization layer for ancestry-probability, hard-boundary, and ancestry-axis map outputs, including a manuscript-style preset inspired by Massatti and Winkler (2022).
