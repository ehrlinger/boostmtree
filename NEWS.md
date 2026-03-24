# boostmtree 2.0.0.9000

## Breaking Changes

* Minimum R version bumped to 4.1.0.
* Maintainership transferred to John Ehrlinger. Original authors (Ishwaran, Pande, Kogalur) remain as authors.

## New Features

* Added `ggRandomForests`, `hvtiPlotR`, and `hvtiRutilities` as optional suggested dependencies for enhanced ggplot2-based visualization workflows.
* Added two vignettes: *Introduction to boostmtree* and *Longitudinal Data Analysis with boostmtree*.
* pkgdown documentation site fully restructured with reference index, articles, and sidebar badges.
* `BugReports` field added to DESCRIPTION pointing to GitHub Issues.

## Documentation

* Comprehensive README rewritten with full badge set, installation instructions, quick-start examples, and reference citations.
* `NEWS.md` introduced (in addition to `inst/NEWS`) for pkgdown changelog rendering.

---

# boostmtree 1.5.2

* Minor maintenance release: updated CI/CD workflows and pkgdown configuration.

# boostmtree 1.5.1

* Added the ability to model a nominal or ordinal longitudinal response. Package now handles continuous, binary, nominal, and ordinal longitudinal responses.

# boostmtree 1.4.1

* Added the ability to model a binary response.
* Added `vimp.boostmtree()` for variable importance estimation using out-of-bag or test data, including joint VIMP for multiple covariates.

# boostmtree 1.3.0

* Fixed `mclapply` on Windows and Debian — defaults to `lapply` on those platforms.
* In-sample cross-validation now uses only out-of-bag data.
* Variable importance (VIMP) calculated from OOB data in grow mode.
* Renamed `vimp.plot` to `vimpPlot`.
* Added `marginalPlot()` — faster alternative to `partialPlot()` providing unadjusted predicted values.

# boostmtree 1.2.1

* Fix to avoid `mclapply` on Debian, defaulting to `lapply` instead.

# boostmtree 1.2.0

* Synchronization with `randomForestSRC` (>= 1.5.0) due to changes in membership output protocols.

# boostmtree 1.1.0

* Added in-sample cross-validation for determining the optimal boosting iteration.
* Various bug fixes.

# boostmtree 1.0.0

* Initial CRAN release.
