# popmaps2

`popmaps2` is the maintained successor to POPMAPS: **Population Management using Ancestry Probability Surfaces**. It estimates spatially explicit ancestry coefficients and ancestry probability surfaces from empirical genetic data across a user-defined landscape.

The package is currently in private alpha. The initial codebase is seeded from the USGS POPMAPS 1.03 release so that results can be compared against the published implementation while the package is modernized, tested, documented, and optimized.

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

The modeling code is known to be computationally expensive. The original implementation loops over raster cells and repeatedly recalculates distances, which can make larger analyses slow. `popmaps2` now includes a faster geographic-distance path for `surface = "G"` that preserves the POPMAPS 1.03 output on validation cases. Least-cost modeling with `surface = "C"` still uses the legacy implementation and remains a priority for modernization.

## Installation

During private development, install from GitHub after authenticating with access to the repository:

```r
install.packages("remotes")
remotes::install_github("landscollective/popmaps2")
```

For local development:

```r
install.packages(c("raster", "sp", "foreach", "doParallel", "gtools", "maps", "plotrix", "MASS", "testthat"))
remotes::install_local(".")
```

Optional legacy functionality currently requires additional packages:

- `gdistance` for least-cost distance surfaces with `surface = "C"`;
- `gplots` for `jackknife_viz()`;
- `viridis` for legacy plotting functions.

`rgeos` has been retired from CRAN and is no longer a package dependency. The `popmap_viz()` boundary plotting path no longer applies the old buffered-boundary adjustment, and the full plotting stack remains a priority target for `terra`/`sf` replacement.

## Data Requirements

The core functions expect three kinds of inputs.

### 1. Raster Surface

`input_raster` defines the extent and resolution of the interpolation. It may be a `terra::SpatRaster`, a legacy `raster::RasterLayer`, or a path readable by `terra::rast()`.

Internally, new input handling is `terra`-first. Some legacy modeling and plotting internals still convert to `raster` objects until those paths are fully modernized.

Raster values are used differently depending on `surface`:

- `surface = "G"` uses geographic distance between empirical sites and raster cells.
- `surface = "C"` uses least-cost distance across raster cell values.

Use `threshold` to skip cells where ancestry coefficients should not be estimated, such as cells below a species distribution model suitability threshold.

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
  threshold = 0,
  ncore = 2
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
| `tune_popmaps()` | Tune geographic-distance POPMAPS parameters with leave-one-out or spatial-block validation metrics. |
| `suggest_tuning_grid()` | Suggest tuning grids from empirical sampling-site distances. |
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

Planned improvements:

- extend `tune_popmaps()` to least-cost surfaces after the `surface = "C"` engine is modernized;
- replace `gdistance` least-cost routines with a maintained alternative;
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
