# popmaps2 0.0.0.9000

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
