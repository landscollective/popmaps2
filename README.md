# popmaps2

`popmaps2` is the maintained successor to POPMAPS: **Population Management
using Ancestry Probability Surfaces**. It estimates spatially explicit ancestry
coefficients and ancestry probability surfaces from empirical ancestry estimates
and user-defined geospatial surfaces.

The package is in development alpha. It is installable, tested, and usable for
current validation work, but it is not on CRAN yet and the public release API may
still change before the first tagged release.

## Scientific Goal

Many restoration and conservation decisions need spatial guidance about
intraspecific genetic variation. `popmaps2` turns empirical ancestry estimates
into raster surfaces that can show:

- hard population boundaries, where each raster cell is assigned to the most
  likely ancestry group;
- ancestry probability, which summarizes confidence in those assignments;
- estimated ancestry coefficients for each ancestry axis or cluster;
- uncertainty where empirical data do not support confident spatial
  interpolation.

The original POPMAPS workflow is described in:

Massatti, R. & Winkler, D. E. (2022). Spatially explicit management of genetic
diversity using ancestry probability surfaces. *Methods in Ecology and
Evolution*, 13, 2668-2681. <https://doi.org/10.1111/2041-210X.13902>

## Current Capabilities

`popmaps2` currently supports:

- the original POPMAPS ancestry-surface workflow;
- faster geographic interpolation for `surface = "G"`;
- internal least-cost distances for `surface = "C"`, replacing the main
  dependency on `gdistance`;
- raster-cell-center interpolation by default, with explicit POPMAPS 1.03
  compatibility mode for historical validation;
- suitability, conductance, and resistance inputs for `surface = "C"`;
- tuning with leave-one-out and spatial-block validation;
- repeated spatial-block validation for uncertainty in validation design;
- candidate-surface comparison with matched validation folds;
- import helpers for raster stacks, point-derived geographic templates, `sf`
  locations, and raster-like EEMS/FEEMS exports;
- report helpers for tuning and surface-comparison diagnostics;
- local validation scripts for empirical example data kept outside the package.

## Scope

`popmaps2` is a downstream ancestry-interpolation and validation tool. It does
not infer population structure, fit species distribution models, run EEMS/FEEMS,
optimize resistance surfaces, or run circuit-theory models.

Those tools are expected to run upstream. Their outputs can then be supplied to
`popmaps2` as empirical ancestry tables or candidate spatial surfaces.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("landscollective/popmaps2")
```

To install source vignettes so `vignette("surface-comparison",
package = "popmaps2")` works from the installed package, build vignettes during
installation:

```r
remotes::install_github(
  "landscollective/popmaps2",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

For local development:

```r
install.packages(c(
  "raster", "sp", "terra", "igraph", "gtools", "maps", "plotrix",
  "MASS", "testthat", "knitr", "rmarkdown", "pkgdown"
))
remotes::install_local(".")
```

Optional legacy functions use optional packages:

- `gdistance` for the legacy `jackknife(surface = "C")` path;
- `gplots` for `jackknife_viz()`;
- `viridis` for legacy plotting palettes.

`rgeos` is retired from CRAN and is not a dependency.

## Input Data

### Raster Surface

`input_raster` defines the interpolation grid. It may be a
`terra::SpatRaster`, a legacy `raster::RasterLayer`, or a file path readable by
`terra::rast()`.

Raster values are interpreted according to `surface`:

| Setting | Meaning |
| --- | --- |
| `surface = "G"` | Geographic interpolation. Raster values are ignored for distances; the raster only supplies geometry and optional prediction masking. |
| `surface = "C"` | Conductance or cost-distance interpolation. Raster values affect movement or gene-flow distance across the landscape. |

For `surface = "C"`, declare the meaning of raster values with
`surface_values`:

| `surface_values` | Meaning |
| --- | --- |
| `"suitability"` | Higher values mean easier movement or stronger support. |
| `"conductance"` | Higher values mean easier movement. |
| `"resistance"` | Higher values mean harder movement; values are inverted internally to conductance. |

`threshold` is a prediction mask. It skips cells where ancestry should not be
estimated, but it is not treated as a movement barrier.

By default, `popmaps()` estimates every output cell at the actual raster cell
center returned by the raster geometry. Set `legacy_compat = TRUE` only when you
need to reproduce POPMAPS 1.03 output exactly for historical comparison. That
legacy mode preserves an old raster-indexing workaround and should not be used
for new analyses.

### Empirical Ancestry Locations

`input_locs` must contain:

| Column | Meaning |
| --- | --- |
| 1 | Sampling location name |
| 2 | Longitude or x-coordinate |
| 3 | Latitude or y-coordinate |
| 4...n | Ancestry coefficients for each ancestry axis or cluster |

Column names may be descriptive, such as `site`, `lon`, `lat`, `axis1`,
`axis2`, and `axis3`, as long as the column order is correct.

Point features can be converted from `sf`:

```r
input_locs <- locs_from_sf(
  sf_points,
  site_col = "site",
  ancestry_cols = c("axis1", "axis2", "axis3")
)
```

### Candidate Surface Helpers

Prepare candidate surfaces explicitly:

```r
geographic <- prepare_popmaps_surface(hija_raster, surface = "G")

sdm_suitability <- prepare_popmaps_surface(
  input_raster = hija_raster,
  surface = "C",
  surface_values = "suitability"
)

resistance <- prepare_popmaps_surface(
  input_raster = hija_raster,
  surface = "C",
  surface_values = "resistance"
)
```

Use converters when inputs are not already single candidate rasters:

```r
point_grid <- surface_from_points(hija_struc, resolution = 0.01)

candidate_surfaces <- surfaces_from_raster_stack(
  input_raster = multi_layer_raster,
  surface = "C",
  surface_values = c("suitability", "conductance", "resistance"),
  include_geographic = TRUE
)

eems_surface <- surface_from_eems(eems_table, value_col = "migration")
feems_surface <- surface_from_feems(feems_table, value_col = "w")
```

`surface_from_points()` defaults to `surface = "G"`. It is mainly a convenience
for building a geographic/template grid from empirical coordinates. A constant
`surface = "C"` grid can be created, but biologically meaningful `C` analyses
should normally come from a supplied suitability, conductance, resistance, EEMS,
FEEMS, or other landscape surface.

## Basic Workflow

Load the package and example data:

```r
library(popmaps2)

data(hija_raster)
data(hija_struc)
data(hija_herb)
```

Aggregate the example raster to keep demonstration runs fast:

```r
ex_raster <- raster::aggregate(hija_raster, fact = 16)
```

Tune geographic-distance parameters:

```r
grid <- suggest_tuning_grid(hija_struc)

tuning <- tune_popmaps(
  input_raster = ex_raster,
  input_locs = hija_struc,
  surface = "G",
  empirical_pt_dist = grid$empirical_pt_dist,
  num_sites = grid$num_sites,
  num_tested = grid$num_tested,
  popmod = grid$popmod,
  quiet = TRUE
)

tuning$best
```

Diagnose how strongly the best parameter combination is supported:

```r
diagnostics <- diagnose_tuning(tuning)

diagnostics$overview
diagnostics$parameter_ranges
```

Use spatial-block validation when the management question involves prediction
into unsampled areas:

```r
spatial_tuning <- tune_popmaps(
  input_raster = ex_raster,
  input_locs = hija_struc,
  validation = "spatial_block",
  n_blocks = 4,
  spatial_block_repeats = 5,
  spatial_block_seed = 1,
  empirical_pt_dist = grid$empirical_pt_dist,
  num_sites = grid$num_sites,
  num_tested = grid$num_tested,
  popmod = grid$popmod,
  quiet = TRUE
)

spatial_tuning$best
```

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
  spatial_block_repeats = 5,
  quiet = TRUE
)

surface_comparison$summary
surface_comparison$support
```

Write a surface-comparison report:

```r
report <- write_surface_comparison_report(
  surface_comparison,
  dir = "surface-comparison-report",
  prefix = "hija-surfaces"
)

report$report
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

For `surface = "C"`, `popmaps2` selects candidate empirical sites by least-cost
distance over the supplied surface. In `legacy_compat = TRUE` mode only, it
reproduces the POPMAPS 1.03 shortcut that first narrowed the candidate pool by
geographic distance before ordering candidates by least-cost distance.

Convert output to a raster and write GeoTIFFs:

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

## Output Structure

`popmaps()` currently returns the original POPMAPS list structure:

| Element | Contents |
| --- | --- |
| `[[1]]` | Hard population boundary matrix |
| `[[2]]` | Ancestry probability matrix |
| `[[3]]...[[n]]` | Estimated ancestry coefficient matrices for each ancestry axis |

The ancestry-axis matrices are weighted ancestry estimates following the
published POPMAPS equation. They are not forced to sum to one at every cell
because the distance-decay weights also carry information about confidence and
distance from empirical data.

Use `popmaps_rast()` and `write_popmaps()` for raster conversion and export.

## Local Validation Scripts

The `tools/` scripts are for local validation, benchmarking, and empirical
example summaries. They are intentionally not run by routine package checks.

All scripts load shared resource settings from `tools/popmaps-script-utils.R`.
They use a conservative fraction of detected logical processors and avoid
package-level `doParallel`, `foreach`, or `parallel` dependencies.

Shared resource environment variables:

| Variable | Default | Purpose |
| --- | --- | --- |
| `POPMAPS_THREADS` | auto | Number of logical processors requested for threaded system libraries. Use `all` for every detected processor. |
| `POPMAPS_THREAD_FRACTION` | `0.75` | Fraction of detected processors used when `POPMAPS_THREADS` is unset. |
| `POPMAPS_TERRA_MEMFRAC` | `0.70` | Fraction of memory `terra` may use before writing temporary files. |
| `POPMAPS_TMPDIR` | R session tempdir | Directory for temporary raster files and intermediate outputs. |

Script-specific thread variables such as `POPMAPS_ASLO_THREADS`,
`POPMAPS_SURFACE_THREADS`, `POPMAPS_EXAMPLE_THREADS`, and
`POPMAPS_TUNING_THREADS` override `POPMAPS_THREADS`.

### Bundled Surface Example

```sh
Rscript tools/example-surface-comparison.R
```

Uses bundled `hija_*` data, creates a point-derived geographic surface, compares
it with the bundled SDM suitability raster, and writes a report under
`local_validation/example_surface_comparison`.

### ASLO Validation

```sh
Rscript tools/validate-aslo-local.R \
  /path/to/aslo_avg.asc \
  /path/to/aslo.txt \
  /tmp/popmaps2-aslo-validation
```

Common controls:

| Variable | Default | Purpose |
| --- | --- | --- |
| `POPMAPS_ASLO_AGGREGATE` | `1` | Aggregate the input raster before modeling. |
| `POPMAPS_ASLO_SURFACE` | `G` | Run geographic or conductance/cost interpolation. |
| `POPMAPS_ASLO_NUM_SITES` | `15` | Set `num_sites`. |
| `POPMAPS_ASLO_NUM_TESTED` | `4` | Set `num_tested`. |
| `POPMAPS_ASLO_POPMOD` | `-0.05` | Set `popmod`. |
| `POPMAPS_ASLO_WRITE_RASTERS` | `false` | Write GeoTIFF output layers. |
| `POPMAPS_ASLO_SAVE_RDS` | `false` | Save the full result object for debugging. |

### Empirical Tuning Examples

```sh
Rscript tools/validate-example-tuning.R ../popmaps_test_data
Rscript tools/summarize-example-tuning.R ../popmaps_test_data/tuning_outputs
```

The first script runs tuning across local `*_avg.asc` and matching `*.txt`
files. The second script summarizes the latest tuning outputs into CSVs,
figures, and `empirical-tuning-report.md`.

### Empirical Surface Comparison

```sh
POPMAPS_SURFACE_AGGREGATE=8 \
POPMAPS_SURFACE_VALIDATION=spatial_block \
POPMAPS_SURFACE_BLOCK_REPEATS=2 \
Rscript tools/compare-example-surfaces.R ../popmaps_test_data
```

This compares a geographic `G` surface against an SDM suitability `C` surface
for each local empirical example. Set
`POPMAPS_SURFACE_INCLUDE_INVERSE=true` to also treat the SDM raster as
resistance for a diagnostic inverse-surface comparison.

### ASLO Benchmarking

```sh
Rscript tools/benchmark-aslo-local.R \
  /path/to/aslo_avg.asc \
  /path/to/aslo.txt \
  /tmp/popmaps2-aslo-benchmark
```

Set `POPMAPS_BENCH_AGGREGATES=16,4,1` and
`POPMAPS_BENCH_SURFACES=G,C` to choose raster sizes and interpolation modes.

## Documentation

Workflow vignettes are available after installing with `build_vignettes = TRUE`:

```r
vignette("parameter-tuning", package = "popmaps2")
vignette("surface-comparison", package = "popmaps2")
vignette("local-empirical-validation", package = "popmaps2")
```

In a repository checkout, the source files are:

```text
vignettes/parameter-tuning.Rmd
vignettes/surface-comparison.Rmd
vignettes/local-empirical-validation.Rmd
```

The repository includes `pkgdown` configuration. To build the website locally:

```r
pkgdown::build_site()
```

The intended GitHub Pages URL is:

<https://landscollective.github.io/popmaps2/>

The pkgdown GitHub Actions workflow is manual-only so the documentation site can
be published without spending Actions minutes on every commit.

## GitHub Actions

The repository uses two workflows:

| Workflow | Trigger | Purpose |
| --- | --- | --- |
| `R-CMD-check` | Pull requests to `main`, plus manual dispatch | Run Linux package checks. macOS checks are manual-only. |
| `pkgdown` | Manual dispatch only | Build and deploy the documentation site to GitHub Pages. |

Routine `R-CMD-check` does not run again on the post-merge push to `main`.
Pull-request checks are the default quality gate, and workflow concurrency
cancels older runs on the same PR or branch.

## Exported Functions

| Function | Purpose |
| --- | --- |
| `popmaps()` | Estimate hard boundaries, ancestry probabilities, and ancestry coefficients across a raster surface. |
| `tune_popmaps()` | Tune geographic or least-cost POPMAPS parameters with leave-one-out or spatial-block validation. |
| `compare_popmaps_surfaces()` | Compare candidate geographic, suitability, conductance, or resistance surfaces with matched validation. |
| `plot_surface_comparison()` | Plot surface validation scores, support gaps, and selected best parameters. |
| `write_surface_comparison_report()` | Write surface-comparison CSVs, diagnostic figures, and a Markdown report. |
| `prepare_popmaps_surface()` | Declare candidate raster semantics before modeling. |
| `surface_from_points()` | Build a simple geographic/template prediction surface from empirical coordinates. |
| `locs_from_sf()` | Convert `sf` point features to a POPMAPS location table. |
| `surfaces_from_raster_stack()` | Convert raster layers to a named candidate-surface list. |
| `surface_from_eems()` | Convert raster-like or coordinate/value EEMS exports to a conductance surface. |
| `surface_from_feems()` | Convert raster-like or coordinate/value FEEMS exports to a conductance surface. |
| `diagnose_tuning()` | Summarize tuning strength, near-best support, and parameter effects. |
| `suggest_tuning_grid()` | Suggest tuning grids from empirical sampling-site distances. |
| `suggest_surface_tuning_grid()` | Suggest tuning grids from distances over a specific candidate surface. |
| `adaptive_tune_popmaps()` | Explore tuning parameter space with random or Latin hypercube sampling and local refinement. |
| `popmaps_rast()` | Convert `popmaps()` list output to a named `terra::SpatRaster`. |
| `write_popmaps()` | Write hard boundary, ancestry probability, and ancestry-axis layers as GeoTIFFs. |
| `anc_extract()` | Extract estimated ancestry coefficients at a coordinate. |
| `jackknife()` | Legacy leave-one-out parameter testing. |
| `jackknife_viz()` | Legacy jackknife heatmap visualization. |
| `popmap_viz()` | Legacy ancestry-surface visualization. |
| `bg_pop_pts()` | Generate and partition random background points by inferred population. |
| `ptsNpop()` | Assign provided sample points to inferred populations. |
| `popmap_pca()` | Build environmental PCA rasters from environmental layers. |
| `envplot()` | Visualize environmental space by inferred population. |

## Relationship To Related Software

`popmaps2` occupies a downstream niche among spatial population-genetic tools.
Common upstream or complementary tools include:

| Software | Primary purpose | Relationship to `popmaps2` |
| --- | --- | --- |
| POPMAPS 1.03 | Original ancestry probability surface package. | Direct predecessor and validation reference. |
| conStruct, LEA, TESS3, ADMIXTURE-style tools | Infer ancestry or spatial population structure. | Potential sources of empirical ancestry estimates. |
| EEMS, FEEMS, reems | Estimate migration or effective-resistance surfaces. | Potential sources of candidate conductance surfaces. |
| ResistanceGA, Circuitscape, Omniscape | Optimize or evaluate resistance/connectivity surfaces. | Potential upstream sources of candidate surfaces or distances. |
| mapmixture, pophelper | Visualize admixture outputs. | Complementary visualization tools, not ancestry interpolation engines. |

## Future Ideas

These are not current package promises. They are candidate directions that should
only move forward as small, testable features:

- precomputed site-site and site-cell distance inputs;
- additional landscape-distance models beyond least-cost distance;
- directional or asymmetric movement surfaces;
- diagnostic plots showing why candidate surfaces differ;
- import helpers for common ancestry-output formats;
- categorical land-cover reclassification into suitability, conductance, or
  resistance values.

## Development Priorities

Near-term priorities are:

1. run empirical validation across all example species with the current
   `G`/`C` surface workflow;
2. benchmark and cache expensive distance calculations;
3. improve modern plotting and map outputs;
4. publish pkgdown documentation;
5. prepare the first tagged development release.

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

The initial `popmaps2` codebase is derived from POPMAPS 1.03, released by the
U.S. Geological Survey. The upstream repository identifies the software as
public domain / CC0-1.0. This repository currently retains CC0 licensing; any
license change would require explicit review.

## Contact

Lands Collective  
<https://www.landscollective.org>  
<rob@landscollective.org>
