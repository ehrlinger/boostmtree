# Introduction to boostmtree

## Overview

**boostmtree** extends Friedman’s (2001) gradient descent boosting
framework to handle *longitudinal* response data by using multivariate
tree base learners and penalized B-splines (P-splines) to model
covariate–time interactions. The package supports four response
families:

| Family       | Description                     |
|:-------------|:--------------------------------|
| `Continuous` | Real-valued outcome (default)   |
| `Binary`     | 0/1 or two-level factor         |
| `Nominal`    | Unordered multi-category factor |
| `Ordinal`    | Ordered multi-category factor   |

Although the package is designed for longitudinal data, it degrades
gracefully to a standard boosting model when each subject has only a
single observation.

## Installation

``` r

# Released version
install.packages("boostmtree")

# Development version
remotes::install_github("ehrlinger/boostmtree")
```

## Simulating longitudinal data

[`simLong()`](https://ehrlinger.github.io/boostmtree/reference/simLong.md)
generates synthetic longitudinal datasets with varying degrees of
covariate–time interaction complexity. Four simulation models are
provided:

| `model` | True formula                                               |
|:-------:|:-----------------------------------------------------------|
|    0    | `y ~ x1 + x3 + x4` (main effects only)                     |
|    1    | `y ~ x1 + x3 + x4 + time * x2` (linear interaction)        |
|    2    | `y ~ x1 + x3 + x4 + time^2 * x2^2` (quadratic interaction) |
|    3    | `y ~ x1 + x3 + exp(x4) + time^2 * x2^2 * x3` (complex)     |

``` r

library(boostmtree)

set.seed(42)
dta <- simLong(
  n      = 50,     # training subjects
  ntest  = 25,     # test subjects
  N      = 5,      # average ~5 time points per subject
  model  = 1,      # linear time × x2 interaction
  family = "Continuous",
  q      = 3       # 3 noise variables added
)

## Data is returned in two formats
## dta$dta   — flat data.frame
## dta$dtaL  — list with $features, $time, $id, $y

cat("Training rows:", sum(dta$dtaL$id %in% unique(dta$dtaL$id)[dta$trn]), "\n")
#> Training rows: 682
cat("True model:   ", dta$f.true, "\n")
#> True model:    y ~ g( x1 + x3 + x4 + I(time * x2) )
cat("Covariates:   ", ncol(dta$dtaL$features), "\n")
#> Covariates:    7
```

## Fitting a model

The main entry point is
[`boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/boostmtree.md).
Key tuning parameters are:

| Parameter | Description                              | Default |
|:----------|:-----------------------------------------|:-------:|
| `M`       | Maximum boosting iterations              |   200   |
| `nu`      | Shrinkage / learning rate                |  0.05   |
| `K`       | Terminal nodes per tree                  |    5    |
| `cv.flag` | Use OOB cross-validation for optimal `M` | `FALSE` |
| `nknots`  | Number of B-spline knots                 |   10    |

``` r

## Extract training indices
trn <- dta$trn

fit <- boostmtree(
  x       = dta$dtaL$features[trn, ],
  tm      = dta$dtaL$time[trn],
  id      = dta$dtaL$id[trn],
  y       = dta$dtaL$y[trn],
  family  = "Continuous",
  M       = 50,
  nu      = 0.05,
  cv.flag = TRUE,
  verbose = FALSE
)

print(fit)
#> boostmtree summary
#> model                       : mtree-Pspline learner 
#> fitting mode                : grow 
#> Family                      : Continuous 
#> number of K-terminal nodes  : 5 
#> regularization parameter    : 0.05 
#> sample size                 : 50 
#> number of variables         : 7 
#> number of unique time points: 15 
#> avg. number of time points  : 9.66 
#> B-spline dimension          : 14 
#> penalization order          : 3 
#> boosting iterations         : 50 
#> optimized number iterations : 50 
#> optimized rho               : 0.8131 
#> optimized phi               : 1.3573 
#> OOB cv RMSE                 : 0.6489
```

The `cv.flag = TRUE` option uses the out-of-bag samples at each
iteration to select an optimal stopping point `Mopt`, protecting against
overfitting without a separate held-out validation set.

## Diagnostic plots

[`plot.boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/plot.boostmtree.md)
produces a multi-panel diagnostic display:

``` r

plot(fit)
```

Panels include:

- **Fitted vs. observed** values (longitudinal profiles)
- **OOB error rate** vs. boosting iteration with the optimal `M` marked
- **Parameter evolution**: smoothing (`lambda`), variance (`phi`),
  correlation (`rho`)

## Predicting on new data

``` r

pred <- predict(
  fit,
  x  = dta$dtaL$features[-trn, ],
  tm = dta$dtaL$time[-trn],
  id = dta$dtaL$id[-trn],
  y  = dta$dtaL$y[-trn]    # optional — enables test-set RMSE
)

print(pred)
#> boostmtree summary
#> model                       : mtree-Pspline learner 
#> fitting mode                : predict 
#> Family                      : Continuous 
#> sample size                 : 25 
#> number of variables         : 7 
#> number of unique time points: 15 
#> avg. number of time points  : 7.96 
#> optimized number iterations : 21 
#> optimized rho               : 0.8246 
#> optimized phi               : 1.2769 
#> test set RMSE               : 0.5838
```

When outcomes are provided,
[`predict()`](https://rdrr.io/r/stats/predict.html) reports the test-set
RMSE (for Continuous) or misclassification rate (for
Binary/Nominal/Ordinal).

``` r

plot(pred)
```

## Variable importance

[`vimp.boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/vimp.boostmtree.md)
computes permutation-based variable importance. For grow objects
(training data) it uses OOB samples; for predict objects it uses the
test set.

``` r

vimp_obj <- vimp.boostmtree(fit)
```

[`vimpPlot()`](https://ehrlinger.github.io/boostmtree/reference/vimpPlot.md)
displays the results. For longitudinal models, positive bars (above the
x-axis) represent main effects and negative bars (below) represent
time-interaction effects.

``` r

vimpPlot(vimp_obj)
```

## Marginal and partial dependence plots

### Marginal plots (fast)

[`marginalPlot()`](https://ehrlinger.github.io/boostmtree/reference/marginalPlot.md)
bins the covariate into quantile groups and plots the *raw* (unadjusted)
predicted mean at each group. It is fast and useful for initial
exploration.

``` r

marginalPlot(fit, xvar.names = c("x1", "x2"), plot.it = TRUE)
```

### Partial dependence plots (adjusted)

[`partialPlot()`](https://ehrlinger.github.io/boostmtree/reference/partialPlot.md)
marginalises over all other covariates at each evaluation point,
providing a confounder-adjusted relationship. It is slower but more
interpretable in the presence of correlated predictors.

``` r

partialPlot(fit, xvar.names = c("x1", "x2"), npts = 10)
```

## Binary response

The workflow is identical for binary outcomes. Simply set
`family = "Binary"` and supply a 0/1 integer or two-level factor as `y`.

``` r

set.seed(7)
dta_bin <- simLong(n = 50, ntest = 0, N = 4, model = 1,
                   family = "Binary", q = 2)

fit_bin <- boostmtree(
  x       = dta_bin$dtaL$features,
  tm      = dta_bin$dtaL$time,
  id      = dta_bin$dtaL$id,
  y       = dta_bin$dtaL$y,
  family  = "Binary",
  M       = 50,
  cv.flag = TRUE,
  verbose = FALSE
)

print(fit_bin)
#> boostmtree summary
#> model                       : mtree-Pspline learner 
#> fitting mode                : grow 
#> Family                      : Binary 
#> number of K-terminal nodes  : 5 
#> regularization parameter    : 0.05 
#> sample size                 : 50 
#> number of variables         : 6 
#> number of unique time points: 12 
#> avg. number of time points  : 6.94 
#> B-spline dimension          : 14 
#> penalization order          : 3 
#> boosting iterations         : 50 
#> optimized number iterations : 50 
#> optimized rho               : 0.2111 
#> optimized phi               : 0.1084 
#> OOB cv RMSE                 : 0.3463
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
