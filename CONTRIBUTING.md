# Contributing

`popmaps2` is in private alpha while the package is stabilized and optimized.

## Development Priorities

- Preserve scientific behavior from POPMAPS 1.03 until tests prove an optimized implementation is equivalent or intentionally different.
- Keep changes small enough to review.
- Add tests before changing model behavior.
- Document any change that affects output values, dependency requirements, or runtime.

## Local Checks

Run these before opening a pull request:

```r
testthat::test_local()
```

From the shell:

```sh
R CMD check --no-manual --as-cran .
```

Some legacy optional features currently require retired or unavailable packages. Those paths should be tested separately until replacements are implemented.
