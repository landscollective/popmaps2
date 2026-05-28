# popmaps2

`popmaps2` is the maintained successor to POPMAPS: **Population Management
using Ancestry Probability Surfaces**. It estimates spatially explicit ancestry
coefficients and ancestry probability surfaces from empirical ancestry estimates
and user-defined geospatial surfaces.

The package is in public alpha. It is installable, tested, documented, and
usable for validation work, but it is not on CRAN yet and the public API may
still change before the software manuscript and CRAN submission.

Use tagged releases when reproducibility matters. The `main` branch may move
quickly while empirical validation, performance work, and figure workflows are
still being refined.

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
- modern map figures for ancestry probability, hard boundaries, and ancestry
  axes, including a manuscript-style preset;
- local validation scripts for empirical example data kept outside the package.

## Scope

`popmaps2` is a downstream ancestry-interpolation and validation tool. It does
not infer population structure, fit species distribution models, run EEMS/FEEMS,
optimize resistance surfaces, or run circuit-theory models.

Those tools are expected to run upstream. Their outputs can then be supplied to
`popmaps2` as empirical ancestry tables or candidate spatial surfaces.

## Known Limitations

The current public alpha is appropriate for package validation and method
development, but users should keep these limits in mind:

- formal empirical validation across additional species and candidate surfaces
  is still ongoing;
- expensive least-cost calculations can still be slow on large, fine-resolution
  rasters;
- Windows checks are not yet part of the routine GitHub Actions gate;
- `jackknife(surface = "C")` remains a legacy compatibility path that still uses
  optional `gdistance`;
- management decisions should be based on biological interpretation, validation
  results, uncertainty, and local expertise rather than a single automated
  "best" model.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("landscollective/popmaps2")
```

To install source vignettes so `vignette("quick-start",
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

## Start Here

The fastest way to learn the package is to run the bundled complete workflow:

```r
library(popmaps2)

start_here <- system.file("examples", "start-here.R", package = "popmaps2")
source(start_here)
```

That script tunes a small grid, runs `popmaps()`, writes GeoTIFF layers, and
writes a manuscript-style PNG map using the bundled `hija_*` example data.

For a narrated version of the same workflow, install vignettes and open:

```r
vignette("quick-start", package = "popmaps2")
```

## Input Data

`input_raster` defines the interpolation grid. It may be a `terra::SpatRaster`,
a legacy `raster::RasterLayer`, or a file path readable by `terra::rast()`.
Raster values are ignored for `surface = "G"` and used as suitability,
conductance, or resistance values for `surface = "C"`.

`input_locs` is an empirical ancestry table with this column order:

| Column | Meaning |
| --- | --- |
| 1 | Sampling location name |
| 2 | Longitude or x-coordinate |
| 3 | Latitude or y-coordinate |
| 4...n | Ancestry coefficients for each ancestry axis or cluster |

Column names may be descriptive as long as the order is correct. Before
modeling starts, `popmaps2` checks that site names are present and
unique, coordinate and ancestry columns are numeric, ancestry coefficients are
non-negative and probability-like, and empirical coordinates fall inside
non-`NA` cells of the raster.

Point features can be converted from `sf`:

```r
input_locs <- locs_from_sf(
  sf_points,
  site_col = "site",
  ancestry_cols = c("axis1", "axis2", "axis3")
)
```

## Core Workflow

`popmaps()` currently returns the original POPMAPS list structure:

| Element | Contents |
| --- | --- |
| `[[1]]` | Hard population boundary matrix |
| `[[2]]` | Ancestry probability matrix |
| `[[3]]...[[n]]` | Estimated ancestry coefficient matrices for each ancestry axis |

Use `popmaps_rast()` and `write_popmaps()` for GeoTIFF conversion and export.
Use `plot_popmaps()` and `write_popmaps_plot()` for modern map figures.
The older `popmap_viz()` function is retained only for POPMAPS 1.03 plotting
compatibility.

The main workflow steps are:

1. prepare an empirical ancestry table and raster surface;
2. tune a biologically plausible parameter grid;
3. optionally compare candidate `G` and `C` surfaces with matched validation;
4. run `popmaps()` with the selected parameters and surface;
5. export rasters, figures, and validation summaries.

## Documentation

Workflow vignettes are available after installing with `build_vignettes = TRUE`:

```r
vignette("quick-start", package = "popmaps2")
vignette("parameter-tuning", package = "popmaps2")
vignette("surface-comparison", package = "popmaps2")
vignette("local-empirical-validation", package = "popmaps2")
```

In a repository checkout, the source files are:

```text
vignettes/quick-start.Rmd
vignettes/parameter-tuning.Rmd
vignettes/surface-comparison.Rmd
vignettes/local-empirical-validation.Rmd
```

The vignettes cover the details that used to make this README too long:

| Vignette | Best for |
| --- | --- |
| `quick-start` | First successful run, output files, and map products. |
| `parameter-tuning` | Tuning grids, spatial-block validation, adaptive tuning, and diagnostics. |
| `surface-comparison` | Testing geographic, suitability, conductance, resistance, EEMS, and FEEMS surfaces. |
| `local-empirical-validation` | Running larger private/local empirical examples outside package checks. |

The repository includes `pkgdown` configuration. To build the website locally:

```r
pkgdown::build_site()
```

The intended GitHub Pages URL is:

<https://landscollective.github.io/popmaps2/>

The pkgdown GitHub Actions workflow is manual-only so the documentation site can
be published without spending Actions minutes on every commit.

## GitHub Actions

GitHub Actions are intentionally conservative to avoid wasting minutes:

| Workflow | Trigger | Purpose |
| --- | --- | --- |
| `R-CMD-check` | Pull requests to `main`, plus manual dispatch | Run Linux package checks. macOS checks are manual-only. |
| `Full platform R-CMD-check` | Manual dispatch only | Run a release-readiness matrix across Ubuntu release, Ubuntu devel, macOS release, and Windows release for a selected branch, tag, or commit. |
| `pkgdown` | Manual dispatch only | Build and deploy the documentation site to GitHub Pages. |

Routine `R-CMD-check` does not run again on the post-merge push to `main`.
Pull-request checks are the default quality gate, and workflow concurrency
cancels older runs on the same PR or branch.
Run the full-platform workflow before public releases or after changes that may
behave differently across operating systems.

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
2. benchmark and cache expensive distance calculations for larger rasters;
3. extend modern map outputs with optional vector overlays, insets, and
   validation-summary callouts;
4. add Windows and optional macOS release checks before CRAN submission;
5. draft the software manuscript with validation, benchmark, and reproducibility
   notes.

## Citation

If you use `popmaps2`, cite the methods paper and the software version or commit
you used. GitHub also reads `CITATION.cff` for software-citation metadata.

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
