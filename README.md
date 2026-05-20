# popmaps2

`popmaps2` is the maintained successor to POPMAPS: **Population Management using Ancestry Probability Surfaces**. It estimates spatially explicit ancestry coefficients and ancestry probability surfaces from empirical genetic data across a user-defined landscape.

The package is currently in development alpha. The initial codebase is seeded from the USGS POPMAPS 1.03 release so that results can be compared against the published implementation while the package is modernized, tested, documented, and optimized.

## Why This Exists

Many restoration and conservation decisions need spatial guidance about intraspecific genetic variation. POPMAPS was designed to translate empirical ancestry assignments into maps that show:

- **Hard population boundaries**, where each raster cell is assigned to the most likely genetic cluster.
- **Ancestry probabilities**, which express confidence in those assignments across the landscape.
- **Uncertainty**, especially where cells are intermediate between sampling locations or influenced by mixed ancestry.

The original workflow is described in:

Massatti, R. & Winkler, D. E. (2022). Spatially explicit management of genetic diversity using ancestry probability surfaces. *Methods in Ecology and Evolution*, 13, 2668-2681. <https://doi.org/10.1111/2041-210X.13902>

## Current Status

This repository is not yet a polished public release. It is a development branch that currently prioritizes:

- preserving the POPMAPS 1.03 scientific workflow;
- making the project installable as a standard R package;
- documenting expected inputs and outputs;
- reducing dependency fragility;
- adding tests before deeper algorithmic refactoring;
- replacing slow or deprecated spatial code with modern `terra`/`sf`-based workflows.

The modeling code is known to be computationally expensive. The original implementation loops over raster cells and repeatedly recalculates distances, which can make larger analyses slow. `popmaps2` now includes faster geographic-distance and least-cost-distance paths for `popmaps(surface = "G")` and `popmaps(surface = "C")`. The `G` path preserves POPMAPS 1.03 output on validation cases, and the modern `C` path uses an internal least-cost helper for suitability, conductance, and resistance surfaces while keeping the original suitability-as-conductance default.

## Relationship to Related Software

`popmaps2` occupies a specific niche among spatial population-genetic tools. It is not intended to replace software that infers population structure, estimates migration surfaces, optimizes resistance surfaces, or visualizes admixture results. Instead, it is a downstream spatial decision-support tool: it takes empirical ancestry estimates and asks how those estimates should be interpolated across a landscape in a way that is biologically defensible, predictively validated, and honest about uncertainty.

The closest conceptual neighbors include:

| Software | Primary purpose | Relationship to `popmaps2` |
| --- | --- | --- |
| [POPMAPS](https://www.usgs.gov/software/popmaps-r-package-estimate-ancestry-probability-surfaces) | Estimate ancestry probability surfaces from empirical ancestry coefficients and raster surfaces. | Direct predecessor. `popmaps2` preserves the original workflow while adding modern package infrastructure, faster geographic interpolation, stronger validation, and clearer tuning outputs. |
| [conStruct](https://rdrr.io/github/gbradburd/conStruct/) | Model continuous and discrete population genetic structure while accounting for spatial covariance. | Strong upstream population-structure model, but not primarily a rasterized ancestry-surface tool for management planning. |
| [EEMS](https://github.com/dipetkov/eems), [FEEMS](https://github.com/NovembreLab/feems), and [reems](https://cran.r-universe.dev/reems) | Estimate effective migration surfaces and spatial variation in gene flow. | Biologically relevant for interpreting spatial genetic structure, but the output is migration or effective resistance rather than ancestry probability surfaces. |
| [ResistanceGA](https://github.com/wpeterman/ResistanceGA) | Optimize landscape resistance surfaces against genetic distances. | Highly relevant to future `surface = "C"` work. `popmaps2` will compare geographic and biologically informed landscape surfaces for ancestry interpolation rather than optimizing resistance surfaces as the final product. |
| [TESS3/tess3r](https://rdrr.io/github/bcm-uga/TESS3_encho_sen/man/tess3r.html), Geneland, LEA, ADMIXTURE-style tools | Infer ancestry coefficients, clusters, or spatial population structure. | Useful upstream sources of empirical ancestry estimates, but not designed to interpolate those estimates across user-defined management rasters. |
| [mapmixture](https://www.rdocumentation.org/packages/mapmixture/versions/1.2.0) and [pophelper](https://www.royfrancis.com/pophelper/) | Visualize admixture or population-structure results. | Complementary visualization tools, not interpolation or tuning frameworks. |
| [assignPOP](https://cran.r-universe.dev/assignPOP/doc/manual.html) | Population assignment and assignment accuracy with cross-validation. | Shares the validation mindset, but focuses on assigning individuals or populations rather than creating continuous ancestry probability surfaces. |
| [adegenet/sPCA](https://rdrr.io/cran/adegenet/man/spca.html) | Exploratory spatial genetic analysis and spatial principal components. | Useful for detecting spatial genetic structure, but not a direct ancestry-surface interpolation workflow. |

The goal of `popmaps2` is therefore to identify an interpolation model that best reflects the spatial genetic structure of a focal species, given available empirical ancestry data. This includes selecting the surface over which ancestry is interpolated, such as geographic distance (`surface = "G"`) or a suitability-weighted landscape surface (`surface = "C"`), and selecting parameters that control how empirical sampling locations contribute to predictions across space. The preferred model should minimize predictive error while avoiding false precision: if the empirical data do not support confident ancestry estimates in some areas, the resulting surfaces should show that uncertainty rather than hide it.

## Model Selection Goal

Parameter tuning is not meant to find universal defaults. The goal is to ask whether a species' empirical ancestry estimates are better predicted by local, broad, sparse, dense, weakly distance-decayed, or strongly distance-decayed interpolation behavior. That choice should be evaluated with withheld empirical sites, and ultimately with competing geographic (`surface = "G"`) and suitability-weighted (`surface = "C"`) surfaces.

A useful model is one that:

- predicts withheld ancestry estimates better than alternative parameter combinations;
- remains honest about poorly supported regions instead of creating false precision;
- produces biologically interpretable distance-decay scales, such as the distances where site weights decay to 50% or 10%;
- shows whether the best model is sharply supported or whether several parameter combinations perform similarly;
- can later be compared across candidate surfaces so model choice reflects spatial genetic structure, dispersal, gene flow, and habitat-mediated connectivity rather than convenience.

See `EMPIRICAL_TUNING_NOTES.md` for the current local empirical-example
interpretation. See `G_VS_C_DESIGN.md` for the intended meaning of geographic
versus suitability-weighted surface comparisons.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("landscollective/popmaps2")
```

To install the source vignette so `vignette("surface-comparison",
package = "popmaps2")` works from the installed package, request vignette
building explicitly:

```r
remotes::install_github(
  "landscollective/popmaps2",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

The vignette source is also available directly in
`vignettes/surface-comparison.Rmd` for repository checkouts.

For local development:

```r
install.packages(c("raster", "sp", "terra", "igraph", "gtools", "maps", "plotrix", "MASS", "testthat", "knitr", "rmarkdown", "pkgdown"))
remotes::install_local(".")
```

Optional legacy plotting and jackknife functionality currently requires
additional packages:

- `gdistance` for the legacy `jackknife(surface = "C")` path;
- `gplots` for `jackknife_viz()`;
- `viridis` for legacy plotting functions.

`rgeos` has been retired from CRAN and is no longer a package dependency. The `popmap_viz()` boundary plotting path no longer applies the old buffered-boundary adjustment, and the full plotting stack remains a priority target for `terra`/`sf` replacement.

## Data Requirements

The core functions expect three kinds of inputs.

### 1. Raster Surface

`input_raster` defines the extent and resolution of the interpolation. It may be a `terra::SpatRaster`, a legacy `raster::RasterLayer`, or a path readable by `terra::rast()`.

Internally, new input handling is `terra`-first. Some legacy modeling and plotting internals still convert to `raster` objects until those paths are fully modernized.

Raster values are used differently depending on `surface`:

- `surface = "G"` uses geographic distance between empirical sites and raster cells. Raster values do not affect distances or ancestry weights, although `NA` values and `threshold` can still define which cells receive estimates.
- `surface = "C"` uses least-cost distance across raster cell values. In the original workflow these values are MaxEnt species distribution model logistic values, interpreted as habitat suitability or occupancy probability. Higher values act as higher conductance/easier movement; lower values increase effective distance.

In short, `G` is a geographic interpolation surface and `C` is a
conductance/cost-distance interpolation surface. `popmaps2` compares how well
these supplied surfaces predict withheld ancestry estimates; it does not infer
the upstream SDM, EEMS, FEEMS, or resistance model.

Use `threshold` to skip cells where ancestry coefficients should not be estimated, such as cells below a species distribution model suitability threshold. `threshold` is a prediction mask rather than a hard least-cost barrier.

Candidate surfaces can be prepared explicitly with:

```r
sdm_surface <- prepare_popmaps_surface(
  input_raster = hija_raster,
  surface = "C",
  surface_values = "suitability"
)

resistance_surface <- prepare_popmaps_surface(
  input_raster = hija_raster,
  surface = "C",
  surface_values = "resistance"
)
```

Convenience converters are available when inputs are not already in the exact
single-raster shape:

```r
point_surface <- surface_from_points(hija_struc, resolution = 0.01)

candidate_surfaces <- surfaces_from_raster_stack(
  input_raster = multi_layer_raster,
  surface = "C",
  surface_values = c("suitability", "conductance", "resistance"),
  include_geographic = TRUE
)

eems_surface <- surface_from_eems(eems_table, value_col = "migration")
feems_surface <- surface_from_feems(feems_table, value_col = "w")
```

`surface_values = "suitability"` and `"conductance"` use raster values directly.
`surface_values = "resistance"` converts values to conductance with an inverse
transform. EEMS/FEEMS-derived gene-flow surfaces should usually be treated as
conductance-like inputs after they are exported to a supported raster or
coordinate/value table. Running SDMs, EEMS/FEEMS, Circuitscape, or ResistanceGA
remains outside the core scope of `popmaps2`; the package focuses on using
candidate surfaces to interpolate ancestry and compare predictive support.

`popmaps(surface = "C")`, `tune_popmaps(surface = "C")`, and
`compare_popmaps_surfaces()` can use suitability, conductance, or resistance
rasters with the internal least-cost distance helper. For `surface = "C"`,
distances are relative cost-distance units rather than geographic kilometers, so
decay summaries should be interpreted as surface-specific distance scales.

### 2. Empirical Genetic Locations

`input_locs` must be a data frame with this structure:

| Column | Meaning |
| --- | --- |
| 1 | Sampling location name |
| 2 | Longitude in decimal degrees |
| 3 | Latitude in decimal degrees |
| 4...n | Ancestry coefficients for each genetic axis or cluster |

The embedded `hija_struc` dataset is the reference format.

The original POPMAPS implementation expected columns named `V1`, `V2`, `V3`, and so on. `popmaps2` now normalizes valid input tables internally, so descriptive column names such as `site`, `lon`, `lat`, `axis1`, `axis2`, and `axis3` are accepted as long as the column order is correct.

Point features can be converted from `sf` when that optional package is
installed:

```r
input_locs <- locs_from_sf(
  sf_points,
  site_col = "site",
  ancestry_cols = c("axis1", "axis2", "axis3")
)
```

`surface_from_points()` defaults to `surface = "G"` and should usually be read
as "make a geographic/template grid around these coordinates." It can create a
constant `C` surface, but biologically meaningful `C` analyses should normally
come from a supplied suitability, conductance, resistance, EEMS, FEEMS, or other
landscape surface.

### 3. Optional Sample Points

Functions such as `ptsNpop()` and `bg_pop_pts()` can assign herbarium, occurrence, or background points to inferred populations. The embedded `hija_herb` dataset provides a simple example.

## Basic Workflow

Load the package and example data:

```r
library(popmaps2)

data(hija_raster)
data(hija_struc)
data(hija_herb)
```

Aggregate the example raster to make the demonstration fast:

```r
ex_raster <- raster::aggregate(hija_raster, fact = 16)
```

Use `tune_popmaps()` to compare geographic-distance parameter combinations with
cross-validation:

```r
grid <- suggest_tuning_grid(hija_struc)

tuning <- tune_popmaps(
  input_raster = ex_raster,
  input_locs = hija_struc,
  surface = "G",
  empirical_pt_dist = grid$empirical_pt_dist,
  num_sites = grid$num_sites,
  num_tested = grid$num_tested,
  popmod = grid$popmod
)

tuning$best
```

The returned object includes `tuning$results`, with one row per parameter
combination, and `tuning$folds`, with one row per withheld sampling site. Results
include `half_distance_km` and `ten_pct_distance_km`, which translate `popmod`
into the distances where ancestry weights decay to 50% and 10% of their initial
value.

Use `diagnose_tuning()` to summarize whether the best combination is strongly
supported or whether several combinations perform nearly as well:

```r
diagnostics <- diagnose_tuning(tuning)

diagnostics$overview
diagnostics$parameter_ranges
```

For a stricter test of whether parameters predict unsampled regions, use
spatial-block validation:

```r
spatial_tuning <- tune_popmaps(
  input_raster = ex_raster,
  input_locs = hija_struc,
  validation = "spatial_block",
  n_blocks = 4,
  empirical_pt_dist = grid$empirical_pt_dist,
  num_sites = grid$num_sites,
  num_tested = grid$num_tested,
  popmod = grid$popmod
)

spatial_tuning$best
```

To reduce dependence on a single block layout, repeat spatial-block validation
with rotated spatial partitions:

```r
repeated_spatial_tuning <- tune_popmaps(
  input_raster = ex_raster,
  input_locs = hija_struc,
  validation = "spatial_block",
  n_blocks = 4,
  spatial_block_repeats = 5,
  spatial_block_seed = 1,
  empirical_pt_dist = grid$empirical_pt_dist,
  num_sites = grid$num_sites,
  num_tested = grid$num_tested,
  popmod = grid$popmod
)

repeated_spatial_tuning$best
```

Repeated spatial-block results include `n_validation_repeats` and
repeat-level standard deviations such as `rmse_repeat_sd`.

Compare candidate geographic and landscape surfaces with matched validation:

```r
candidate_surfaces <- list(
  geographic = prepare_popmaps_surface(ex_raster, surface = "G"),
  suitability = prepare_popmaps_surface(
    ex_raster,
    surface = "C",
    surface_values = "suitability"
  )
)

surface_comparison <- compare_popmaps_surfaces(
  input_locs = hija_struc,
  surfaces = candidate_surfaces,
  validation = "spatial_block",
  spatial_block_repeats = 5
)

surface_comparison$summary
surface_comparison$support
surface_comparison$grids
```

Draw quick diagnostics or write durable report artifacts:

```r
plot_surface_comparison(surface_comparison, type = "score")
plot_surface_comparison(surface_comparison, type = "percent_from_best")

report <- write_surface_comparison_report(
  surface_comparison,
  dir = "surface-comparison-report",
  prefix = "hija-surfaces"
)

report$report
```

This comparison asks which supplied surface best predicts withheld empirical
ancestry estimates under the POPMAPS interpolation workflow. It does not replace
upstream landscape-genetic or SDM analyses. By default, missing `popmod` and
`empirical_pt_dist` values are suggested separately for each surface from that
surface's empirical site-distance matrix. This avoids forcing least-cost
distances through a tuning grid scaled for geographic kilometers. Set
`surface_grid = "shared"` or pass explicit `popmod` and `empirical_pt_dist`
values when a deliberately shared grid is desired.

For local empirical examples stored outside the package repository, run:

```sh
POPMAPS_SURFACE_AGGREGATE=8 \
POPMAPS_SURFACE_VALIDATION=spatial_block \
POPMAPS_SURFACE_BLOCK_REPEATS=2 \
Rscript tools/compare-example-surfaces.R
```

The script looks for `*_avg.asc` and matching `*.txt` files in
`../popmaps_test_data`, compares geographic (`G`) and SDM suitability (`C`)
surfaces, and writes aggregate summary CSVs plus one standard
`write_surface_comparison_report()` folder per species/validation comparison.

For a minimal bundled-data example using the input converters and report helper:

```sh
Rscript tools/example-surface-comparison.R
```

This creates a geographic surface from the embedded example coordinates with
`surface_from_points()`, compares it to the embedded SDM suitability raster, and
writes a report under `local_validation/example_surface_comparison`.

A longer walkthrough is available in the vignette:

```r
vignette("surface-comparison", package = "popmaps2")
```

For GitHub/source installs, build vignettes during installation or read the
source file at `vignettes/surface-comparison.Rmd`.

The legacy `jackknife()` function is still available for compatibility with
`jackknife_viz()`.

For larger candidate spaces, use `adaptive_tune_popmaps()` to sample parameter
space and refine around the best-performing region:

```r
adaptive <- adaptive_tune_popmaps(
  input_raster = ex_raster,
  input_locs = hija_struc,
  n_initial = 50,
  n_refine = 50,
  seed = 1
)

adaptive$best
```

Estimate an ancestry probability surface:

```r
aps <- popmaps(
  input_raster = ex_raster,
  input_locs = hija_struc,
  surface = "G",
  empirical_pt_dist = 5,
  num_sites = 15,
  num_tested = 4,
  popmod = -0.05,
  threshold = 0
)
```

Extract ancestry coefficients for a location:

```r
anc_extract(
  pop_raster_list = aps,
  input_raster = ex_raster,
  input_locs = hija_struc,
  dec_lat = 39.46522,
  dec_long = -110.9525
)
```

Convert the output to a `terra::SpatRaster` and write GeoTIFF layers:

```r
aps_raster <- popmaps_rast(aps, ex_raster)

write_popmaps(
  pop_raster_list = aps,
  input_raster = ex_raster,
  dir = "outputs",
  prefix = "hija",
  overwrite = TRUE
)
```

Run the built-in POPMAPS 1.03 baseline validation:

```r
validate_popmaps_baseline()
```

This validation uses the bundled `Hilaria jamesii` example data, runs a small ancestry probability surface, and compares every output matrix to a frozen POPMAPS 1.03 reference result. It is intended to catch unintended scientific drift before deeper optimization work.

## Output Structure

`popmaps()` currently returns a list that preserves the original POPMAPS structure:

| Element | Contents |
| --- | --- |
| `[[1]]` | Hard population boundary matrix |
| `[[2]]` | Ancestry probability matrix |
| `[[3]]...[[n]]` | Estimated ancestry coefficient matrices for each genetic axis |

A future release will wrap this list in an S3 class with helper methods for printing, plotting, raster conversion, and export.

## Exported Functions

| Function | Purpose |
| --- | --- |
| `popmaps()` | Estimate hard boundaries, ancestry probabilities, and ancestry coefficients across a raster surface. |
| `tune_popmaps()` | Tune geographic or least-cost POPMAPS parameters with leave-one-out or spatial-block validation metrics. |
| `compare_popmaps_surfaces()` | Compare candidate geographic, suitability, conductance, or resistance surfaces with matched validation. |
| `plot_surface_comparison()` | Plot surface validation scores, relative support gaps, and selected best parameters. |
| `write_surface_comparison_report()` | Write surface-comparison CSVs, diagnostic figures, and a Markdown report. |
| `surface_from_points()` | Build a simple prediction surface from empirical coordinates. |
| `locs_from_sf()` | Convert `sf` point features to a POPMAPS location table. |
| `surfaces_from_raster_stack()` | Convert each raster layer to a candidate surface list. |
| `surface_from_eems()` | Convert raster-like or coordinate/value EEMS exports to a conductance surface. |
| `surface_from_feems()` | Convert raster-like or coordinate/value FEEMS exports to a conductance surface. |
| `diagnose_tuning()` | Summarize tuning strength, near-best support, and parameter effects. |
| `suggest_tuning_grid()` | Suggest tuning grids from empirical sampling-site distances. |
| `suggest_surface_tuning_grid()` | Suggest tuning grids from distances measured over a geographic or least-cost candidate surface. |
| `adaptive_tune_popmaps()` | Explore tuning parameter space with random or Latin hypercube sampling and local refinement. |
| `jackknife()` | Test parameter combinations with a leave-one-out approach. |
| `jackknife_viz()` | Visualize jackknife performance as heatmaps. |
| `anc_extract()` | Extract estimated ancestry coefficients at a geographic coordinate. |
| `popmap_viz()` | Legacy visualization of ancestry probability surfaces. |
| `popmaps_rast()` | Convert `popmaps()` list output to a named `terra::SpatRaster`. |
| `write_popmaps()` | Write hard boundary, ancestry probability, and ancestry-axis layers as GeoTIFFs. |
| `bg_pop_pts()` | Generate and partition random background points by inferred population. |
| `ptsNpop()` | Assign provided sample points to inferred populations. |
| `popmap_pca()` | Build environmental PCA rasters from environmental layers. |
| `envplot()` | Visualize environmental space by inferred population. |

## Documentation Website

The repository now includes a `pkgdown` scaffold for a future public reference
site. To preview it locally after installing development dependencies, run:

```r
pkgdown::build_site()
```

The intended public URL is <https://landscollective.github.io/popmaps2/> unless
the project later moves under a Lands Collective organization account.

## Larger Local Validation

Real project datasets should generally stay outside the package repository unless they are cleared for redistribution. The ASLO validation workflow is therefore a local script that points at files on your machine and writes summaries to a temporary or user-specified output directory.

Example:

```sh
R CMD INSTALL .

Rscript tools/validate-aslo-local.R \
  /path/to/aslo_avg.asc \
  /path/to/aslo.txt \
  /tmp/popmaps2-aslo-validation
```

By default the script writes:

- `aslo-run-summary.csv`, with parameters and runtime;
- `aslo-layer-summary.csv`, with dimensions, non-NA counts, and value ranges by output layer.

Optional environment variables:

| Variable | Default | Purpose |
| --- | --- | --- |
| `POPMAPS_ASLO_AGGREGATE` | `1` | Aggregate the input raster before modeling. |
| `POPMAPS_ASLO_NUM_SITES` | `15` | Set `num_sites`. |
| `POPMAPS_ASLO_NUM_TESTED` | `4` | Set `num_tested`. |
| `POPMAPS_ASLO_POPMOD` | `-0.05` | Set `popmod`. |
| `POPMAPS_ASLO_THRESHOLD` | `0` | Set `threshold`. |
| `POPMAPS_ASLO_WRITE_RASTERS` | `false` | Write GeoTIFF output layers under `rasters/`. |
| `POPMAPS_ASLO_SAVE_RDS` | `false` | Save the full R result object for debugging. |

The local validation scripts also share resource settings that are intended to
work across macOS, Linux, Windows, and common cluster environments without
requiring `doParallel`, `foreach`, or package-level use of R's `parallel`
package:

| Variable | Default | Purpose |
| --- | --- | --- |
| `POPMAPS_THREADS` | auto | Set logical processors used by threaded system libraries. Use `all` to request every detected processor. |
| `POPMAPS_THREAD_FRACTION` | `0.75` | Fraction of detected processors used when `POPMAPS_THREADS` is unset. One processor is left free on machines with more than two processors. |
| `POPMAPS_TERRA_MEMFRAC` | `0.70` | Fraction of memory `terra` may use before writing temporary files. |
| `POPMAPS_TMPDIR` | R session tempdir | Directory for temporary raster files and other temporary outputs. |

Script-specific overrides such as `POPMAPS_ASLO_THREADS`,
`POPMAPS_EXAMPLE_THREADS`, `POPMAPS_SURFACE_THREADS`, and
`POPMAPS_TUNING_THREADS` take precedence over `POPMAPS_THREADS`.

To rerun tuning validation across empirical examples kept outside the package,
place `*_avg.asc` rasters and matching `*.txt` location files in a directory and
run:

```sh
Rscript tools/validate-example-tuning.R ../popmaps_test_data
```

The script writes best-parameter summaries, near-best parameter support, and
parameter-effect tables. By default it runs exhaustive grid tuning for both
leave-one-site-out and spatial-block validation. Set
`POPMAPS_EXAMPLE_SEARCH=adaptive` to use adaptive sampling instead. Set
`POPMAPS_EXAMPLE_BLOCK_REPEATS=5` and `POPMAPS_EXAMPLE_VALIDATION=spatial_block`
to run repeated spatial-block validation.

For runtime and memory-scale checks on ASLO-like inputs, run:

```sh
Rscript tools/benchmark-aslo-local.R /path/to/aslo_avg.asc /path/to/aslo.txt /tmp/popmaps2-aslo-benchmark
```

Set `POPMAPS_BENCH_AGGREGATES=16,4,1` and `POPMAPS_BENCH_SURFACES=G,C` to
control the raster sizes and surfaces included in the benchmark. The output CSV
records elapsed time, R object memory, garbage-collector maximum memory, raster
dimensions, and resource settings.

Summarize the latest empirical tuning run with:

```sh
Rscript tools/summarize-example-tuning.R ../popmaps_test_data/tuning_outputs
```

This creates a timestamped report directory with `empirical-tuning-report.md`,
summary CSVs, and figures for best validation score, near-best support,
distance-decay scales, and parameter effects. When repeated spatial-block output
is available, the best-score plot includes repeat-level error bars.

## Future Ideas Sandbox

This section is a parking place for ideas that may be useful later but are not
current package promises. The core scope remains: take empirical ancestry
coefficients plus user-defined spatial inputs, interpolate ancestry surfaces,
validate parameter choices, and report uncertainty clearly.

Potentially useful future directions:

- **Candidate-surface comparison**: allow users to provide several surfaces and
  compare their predictive support for ancestry interpolation with matched
  validation folds. This would be predictive model selection for POPMAPS
  surfaces, not a replacement for full landscape-genetic hypothesis testing.
- **Precomputed distance inputs**: allow advanced users to supply site-site and
  site-cell distance matrices directly. This would make `popmaps2` compatible
  with distances generated by Circuitscape, Omniscape, EEMS/FEEMS-derived
  workflows, custom graph models, or other future tools without bringing those
  modeling engines into the package.
- **Additional landscape-distance models**: least-cost distance is the legacy
  `surface = "C"` model, but other biological hypotheses may be better for some
  species. Effective-resistance, commute-distance, or randomized-shortest-path
  approaches could represent gene flow through many parallel routes rather than
  only the single easiest route.
- **Directional or asymmetric movement**: some systems may need transitions that
  differ by direction, such as rivers, slopes, wind, ocean currents, or
  directional dispersal. This would require a more general graph interface than
  the current POPMAPS 1.03-compatible conductance model.
- **Surface diagnostics**: plots or summaries could show where candidate
  surfaces create similar or very different distance relationships, which may
  help users interpret why one surface validates better than another.
- **Import helpers for upstream ancestry outputs**: lightweight converters for
  STRUCTURE, ADMIXTURE, LEA, TESS3, conStruct, or similar ancestry outputs could
  make it easier to create the `input_locs` table expected by `popmaps2`.
- **Categorical land-cover reclassification**: users could provide land-cover
  rasters plus an explicit table translating classes into suitability,
  conductance, or resistance values.

Ideas in this sandbox should only move into the roadmap after they can be stated
as small, testable features that support ancestry-surface interpolation rather
than expanding `popmaps2` into a general SDM, landscape-genetics, or circuit
theory package.

## Optimization Plan

The highest-priority performance work is in `popmaps()` and `jackknife()`.

Completed:

- precompute empirical-site distance matrices instead of recalculating distances repeatedly;
- vectorize geographic-distance calculations;
- avoid repeated `raster::extract()` calls inside geographic-distance cell loops;
- test optimized geographic outputs against POPMAPS 1.03 reference outputs.
- add a fast geographic-distance tuning workflow that scores parameter combinations at withheld empirical sites.
- suggest data-adaptive tuning grids from empirical site distances and sample larger parameter spaces adaptively.
- report biologically interpretable distance-decay scales and support spatial-block tuning validation.
- summarize tuning strength, near-best parameter support, and parameter effects with `diagnose_tuning()`.
- add a repeatable local empirical-example tuning validation script and reporting workflow.
- support repeated spatial-block validation with repeat-level uncertainty summaries.
- add `prepare_popmaps_surface()` to declare whether candidate `C` rasters
  represent suitability, conductance, or resistance before modern least-cost
  modeling.
- add internal least-cost distance helpers that match the legacy `gdistance`
  path on small validation rasters.
- extend `tune_popmaps()` to suitability-, conductance-, and
  resistance-weighted least-cost surfaces.
- route full `popmaps(surface = "C")` ancestry-surface estimation through the
  internal least-cost helper.
- add `compare_popmaps_surfaces()` for matched predictive comparison of
  user-supplied geographic and landscape surfaces.

Planned improvements:

- replace `raster`, `sp`, and `rgeos` plotting internals with `terra` and `sf`;
- add progress reporting and reproducible parallel execution;
- benchmark legacy and optimized implementations on small, medium, and full-size rasters;
- extend validation coverage with larger real-world datasets.

## Development Roadmap

### Phase 0: Baseline

- Import POPMAPS 1.03 source.
- Rename package to `popmaps2`.
- Document provenance and current limitations.
- Add local package infrastructure.

### Phase 1: Reliable R Package

- Add tests for data loading, argument validation, and small example outputs.
- Update roxygen documentation.
- Add GitHub Actions for `R CMD check`.
- Add a vignette that reproduces the published `Hilaria jamesii` example.
- Add a `pkgdown` documentation site.

### Phase 2: Performance Modernization

- Profile `popmaps()` and `jackknife()`.
- Implement a tested fast geographic-distance path. `(initial surface = "G" engine complete)`
- Implement a maintained least-cost path.
- Add benchmark results to the documentation.

### Phase 3: Public Release

- Prepare a tagged release.
- Publish package documentation on the Lands Collective website.
- Archive a reproducible release if needed.
- Prepare a software release manuscript.

## Citation

If you use `popmaps2`, cite the methods paper and the software version you used.

```bibtex
@article{Massatti2022AncestryProbabilitySurfaces,
  author = {Massatti, Rob and Winkler, Daniel E.},
  title = {Spatially explicit management of genetic diversity using ancestry probability surfaces},
  journal = {Methods in Ecology and Evolution},
  volume = {13},
  pages = {2668--2681},
  year = {2022},
  doi = {10.1111/2041-210X.13902}
}
```

The original POPMAPS 1.03 release is available from USGS GitLab:

<https://code.usgs.gov/GWRC/popmaps>

## License And Provenance

The initial `popmaps2` codebase is derived from POPMAPS 1.03, released by the U.S. Geological Survey. The upstream repository identifies the software as public domain / CC0-1.0. This repository retains CC0 licensing unless a future release explicitly changes that after review.

## Contact

Lands Collective  
<https://www.landscollective.org>  
<rob@landscollective.org>
