# popmaps2 Roadmap

This file tracks the practical path from the POPMAPS 1.03 baseline to a public `popmaps2` R package and software manuscript.

## Near Term

- Keep the POPMAPS 1.03 implementation available as the reference baseline.
- Add tests around the embedded `Hilaria jamesii` datasets.
- Continue moving public inputs through the new `terra`-first validation layer.
- Update examples so they run quickly and do not require retired packages.
- Add a vignette that reproduces the published workflow at a reduced raster resolution.
- Profile `popmaps()` and `jackknife()` on representative rasters.

## Performance Work

- Precompute empirical-site distances.
- Vectorize geographic-distance calculations.
- Reduce repeated raster extraction inside nested loops.
- Add a single-cell or small-grid internal estimator that can be unit tested.
- Benchmark serial and parallel execution.
- Compare optimized results to the baseline using tolerances documented in tests.

## Spatial Modernization

- Replace `raster` internals with `terra` where practical.
- Replace `sp` objects with `sf`/matrix/data-frame interfaces where practical.
- Replace `rgeos::gBuffer()` in plotting.
- Identify a maintained replacement for `gdistance` least-cost workflows.

## Release Readiness

- Pass `R CMD check` on macOS, Linux, and Windows.
- Build a `pkgdown` site.
- Add lifecycle badges and a changelog.
- Create a tagged GitHub release.
- Publish package documentation on the Lands Collective website.
- Draft software release manuscript with benchmarks and reproducibility notes.
