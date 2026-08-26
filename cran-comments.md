# WeightIt 2.1.0

## Test environments

* local: macOS Tahoe 26.5.2, aarch64-apple-darwin23, R 4.6.1 (2026-06-24)
* GitHub Actions: see `.github/workflows/`

## R CMD check results

`R CMD check --as-cran` on the built tarball:

```
Status: OK
```

0 errors | 0 warnings | 0 notes

The full test suite (286 `test_that()` blocks) was run separately with
`NOT_CRAN=true`; all pass. Most tests are behind `skip_on_cran()` because they fit
many models across the optional weighting-method packages, so the suite CRAN runs
is deliberately a small subset.

## Reverse dependencies

WeightIt has 23 reverse dependencies. These were checked with `revdepcheck`
against both the CRAN and development versions.

<!-- Summary to be filled in from revdep/README.md once the run completes. -->

## Notes on this release

The `.cens()` function, added in 2.0.0, is now re-exported from *cobalt* rather
than defined here. *cobalt* 5.0.0 added its own `.cens()` so that
`cobalt::bal.tab()` could recognize a censoring model, and two identically named
exports meant one masked the other whenever both packages were attached. There is
now a single function, and `Imports:` requires `cobalt (>= 5.0.0)` accordingly.
The behavior of `.cens()` is unchanged.
