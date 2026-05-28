# boostmtree — Boosted Multivariate Trees for Longitudinal Data

**boostmtree** implements Friedman’s (2001) gradient descent boosting
algorithm for modeling longitudinal responses using multivariate tree
base learners. Covariate-time interactions are captured via penalized
B-splines (P-splines) with an adaptively estimated smoothing parameter.
The package handles continuous, binary, nominal, and ordinal responses,
and works equally well for cross-sectional data.

## Features

- Gradient boosting with multivariate tree base learners (powered by
  [randomForestSRC](https://cran.r-project.org/package=randomForestSRC))
- P-spline modeling of covariate × time interactions with adaptive
  smoothing
- Stochastic gradient descent with bootstrap sampling and OOB error
  estimation
- In-sample cross-validation for optimal stopping iteration selection
- Variable importance (VIMP) via permutation, including joint VIMP
- Partial and marginal dependence plots
- Supports **Continuous**, **Binary**, **Nominal**, and **Ordinal**
  longitudinal responses
- Parallel processing via the `parallel` package
- Missing covariate imputation (on-the-fly, via `randomForestSRC`)

## Installation

Install the released version from CRAN:

``` r

install.packages("boostmtree")
```

Install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("ehrlinger/boostmtree")
```

### Optional enhanced visualization

For ggplot2-based visualization workflows, install the suggested
packages:

``` r

remotes::install_github("ehrlinger/ggRandomForests")
remotes::install_github("ehrlinger/hvtiPlotR")
remotes::install_github("ehrlinger/hvtiRutilities")
```

## Quick Start

### Simulated continuous longitudinal data

``` r

library(boostmtree)

## Simulate training and test data (model 1: linear time × covariate interaction)
set.seed(42)
dta <- simLong(n = 100, ntest = 100, model = 1, family = "Continuous", q = 5)

## Fit the model
fit <- boostmtree(
  x  = dta$dtaL$features[dta$trn, ],
  tm = dta$dtaL$time[dta$trn],
  id = dta$dtaL$id[dta$trn],
  y  = dta$dtaL$y[dta$trn],
  family = "Continuous",
  M  = 200,
  cv.flag = TRUE   # enables OOB-based optimal stopping
)

print(fit)
plot(fit)
```

### Predict on held-out data

``` r

pred <- predict(
  fit,
  x  = dta$dtaL$features[-dta$trn, ],
  tm = dta$dtaL$time[-dta$trn],
  id = dta$dtaL$id[-dta$trn],
  y  = dta$dtaL$y[-dta$trn]
)

print(pred)
plot(pred)
```

### Variable importance

``` r

## Individual VIMP for all covariates
vimp_obj <- vimp.boostmtree(fit)
vimpPlot(vimp_obj)

## Joint VIMP for a pair of variables
vimp_joint <- vimp.boostmtree(fit, x.names = c("x1", "x2"), joint = TRUE)
```

### Partial and marginal dependence plots

``` r

## Marginal plot (fast, unadjusted)
marginalPlot(fit, xvar.names = c("x1", "x2"))

## Partial dependence plot (slower, confounder-adjusted)
partialPlot(fit, xvar.names = "x1")
```

### Binary response — Atrial Fibrillation data

``` r

data(AF, package = "boostmtree")

fit_af <- boostmtree(
  x      = AF$feature,
  tm     = AF$time,
  id     = AF$id,
  y      = AF$y,
  family = "Binary",
  M      = 300,
  cv.flag = TRUE
)

print(fit_af)
plot(fit_af)

## Variable importance
vimp_af <- vimp.boostmtree(fit_af)
vimpPlot(vimp_af)
```

### Continuous response — Spirometry data

``` r

data(spirometry, package = "boostmtree")

fit_spi <- boostmtree(
  x      = spirometry$features,
  tm     = spirometry$time,
  id     = spirometry$id,
  y      = spirometry$y,
  family = "Continuous",
  M      = 300,
  cv.flag = TRUE
)

print(fit_spi)
marginalPlot(fit_spi, xvar.names = colnames(spirometry$features)[1:4])
```

## Key Functions

| Function | Description |
|----|----|
| [`boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/boostmtree.md) | Fit a boosted multivariate tree model |
| [`predict.boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/predict.boostmtree.md) | Predict on new data |
| [`print.boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/print.boostmtree.md) | Print model summary |
| [`plot.boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/plot.boostmtree.md) | Diagnostic plots (error curve, fitted values, residuals) |
| [`vimp.boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/vimp.boostmtree.md) | Variable importance scores |
| [`vimpPlot()`](https://ehrlinger.github.io/boostmtree/reference/vimpPlot.md) | Plot variable importance |
| [`marginalPlot()`](https://ehrlinger.github.io/boostmtree/reference/marginalPlot.md) | Marginal dependence plots (fast) |
| [`partialPlot()`](https://ehrlinger.github.io/boostmtree/reference/partialPlot.md) | Partial dependence plots (adjusted) |
| [`simLong()`](https://ehrlinger.github.io/boostmtree/reference/simLong.md) | Simulate longitudinal data |

## Datasets

| Dataset | Subjects | Observations | Response | Description |
|----|----|----|----|----|
| `AF` | 228 | 7,949 | Binary | Atrial fibrillation presence/absence after surgical ablation |
| `spirometry` | 509 | 9,471 | Continuous | FEV1% after lung transplantation |

## Documentation

Full documentation is available at
**<https://ehrlinger.github.io/boostmtree/>**, including:

- [Reference manual](https://ehrlinger.github.io/boostmtree/reference/)
- [Introduction
  vignette](https://ehrlinger.github.io/boostmtree/articles/introduction.html)
- [Longitudinal analysis
  vignette](https://ehrlinger.github.io/boostmtree/articles/longitudinal-analysis.html)

## Authors

- **Hemant Ishwaran** — original methodology
- **Amol Pande** — original methodology
- **Udaya B. Kogalur** — original implementation
- **John Ehrlinger** — current maintainer

## Citation

If you use **boostmtree** in your research, please cite:

> Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
> E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
> data. *Machine Learning*, 106(2): 277–305. doi:
> [10.1007/s10994-016-5597-1](https://doi.org/10.1007/s10994-016-5597-1)

``` r

citation("boostmtree")
```

## References

Friedman J.H. (2001). Greedy function approximation: a gradient boosting
machine. *Annals of Statistics*, 29(5): 1189–1232.

Friedman J.H. (2002). Stochastic gradient boosting. *Computational
Statistics & Data Analysis*, 38(4): 367–378.

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
data. *Machine Learning*, 106(2): 277–305. doi:
[10.1007/s10994-016-5597-1](https://doi.org/10.1007/s10994-016-5597-1)

## License

GPL (\>= 3) — see
[LICENSE](https://ehrlinger.github.io/boostmtree/LICENSE) for details.
