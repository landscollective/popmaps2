# popmaps2 Roadmap

This file tracks the practical path from the POPMAPS 1.03 baseline to a public `popmaps2` R package and software manuscript.

## Near Term

- Keep the POPMAPS 1.03 implementation available as the reference baseline.
- Add tests around the embedded `Hilaria jamesii` datasets.
- Continue moving public inputs through the new `terra`-first validation layer.
- Add user-facing raster conversion and GeoTIFF export helpers.
- Add a modern parameter-tuning helper for geographic-distance leave-one-site-out validation.
- Add data-adaptive tuning-grid suggestions and adaptive parameter-space search.
- Add spatial-block validation and biologically interpretable distance-decay summaries for tuning.
- Add tuning diagnostics that distinguish strong parameter support from broad near-best support.
- Add repeated spatial-block validation so parameter support can be evaluated across multiple spatial partitions.
- Document the intended `surface = "G"` versus `surface = "C"` model-selection contract before modernizing least-cost distance code.
- Add a surface-preparation object that records whether candidate rasters are suitability, conductance, or resistance inputs.
- Maintain local validation scripts for larger datasets that should not be committed.
- Keep empirical example validation repeatable from local, uncommitted `*_avg.asc` and `*.txt` files.
- Update examples so they run quickly and do not require retired packages.
- Add a vignette that reproduces the published workflow at a reduced raster resolution.
- Profile `popmaps()` and `jackknife()` on representative rasters.

## Performance Work

- Precompute empirical-site distances.
- Vectorize geographic-distance calculations.
- Reduce repeated raster extraction inside nested loops.
- Add a single-cell or small-grid internal estimator that can be unit tested.
- Add internal least-cost distance helpers and validate them against the legacy `gdistance` path on small rasters.
- Route full `popmaps(surface = "C")` interpolation through the internal least-cost helper while preserving POPMAPS 1.03 suitability-as-conductance behavior.
- Extend parameter tuning to modernized suitability-, conductance-, and resistance-weighted least-cost surfaces.
- Evaluate `surface = "G"` and `surface = "C"` with matched validation folds, metrics, and uncertainty summaries.
- Add candidate-surface comparison for user-supplied rasters without running upstream SDM, EEMS/FEEMS, Circuitscape, or ResistanceGA models inside `popmaps2`. `(surface-specific grid and empirical-report workflow complete)`
- Run sensitivity checks for empirical surface comparison across raster aggregation, spatial-block repeat count, and near-best tolerance.
- Benchmark serial and parallel execution.
- Compare optimized results to the baseline using tolerances documented in tests.

## Spatial Modernization

- Replace `raster` internals with `terra` where practical.
- Replace `sp` objects with `sf`/matrix/data-frame interfaces where practical.
- Replace `rgeos::gBuffer()` in plotting.
- Finish replacing `gdistance` in remaining legacy compatibility paths, especially `jackknife(surface = "C")`.

## Release Readiness

- Pass `R CMD check` on macOS, Linux, and Windows.
- Build a `pkgdown` site.
- Add lifecycle badges and a changelog.
- Create a tagged GitHub release.
- Publish package documentation on the Lands Collective website.
- Draft software release manuscript with benchmarks and reproducibility notes.
