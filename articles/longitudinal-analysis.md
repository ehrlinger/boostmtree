# Longitudinal Data Analysis with boostmtree

## Motivation

Standard regression methods struggle with longitudinal data: repeated
measures are correlated within subjects, observation schedules are often
unbalanced, and covariate effects can change over time. **boostmtree**
handles all three by combining multivariate tree base learners
(non-linear covariate effects) with adaptive P-splines (covariate–time
interactions), without requiring imputation or a balanced panel.

This vignette demonstrates the full analysis workflow using the two
bundled real datasets.

------------------------------------------------------------------------

## Spirometry: Continuous longitudinal response

The `spirometry` dataset contains 9,471 post-lung-transplant FEV1%
measurements from 509 patients along with 23 baseline clinical
covariates.

``` r

library(boostmtree)

data(spirometry, package = "boostmtree")

cat("Patients:         ", length(unique(spirometry$id)), "\n")
#> Patients:          509
cat("Total obs:        ", length(spirometry$y), "\n")
#> Total obs:         9471
cat("Covariates:       ", ncol(spirometry$features), "\n")
#> Covariates:        23
cat("FEV1% range:      ", round(range(spirometry$y), 1), "\n")
#> FEV1% range:       9.5 165.3
cat("Time range (yrs): ", round(range(spirometry$time), 2), "\n")
#> Time range (yrs):  0.01 15.42
```

### Data structure

[`boostmtree()`](https://ehrlinger.github.io/boostmtree/reference/boostmtree.md)
expects data in *long format* (one row per observation), split into four
components:

| Argument | Description                                               |
|:---------|:----------------------------------------------------------|
| `x`      | `n_obs × p` covariate matrix (time-invariant per subject) |
| `tm`     | Numeric vector of observation times (length `n_obs`)      |
| `id`     | Subject identifier vector (length `n_obs`)                |
| `y`      | Response vector (length `n_obs`)                          |

``` r

set.seed(2024)

## Train/test split by subject
subjects  <- unique(spirometry$id)
n_train   <- floor(0.8 * length(subjects))
trn_subj  <- sample(subjects, n_train)
trn_idx   <- spirometry$id %in% trn_subj

cat("Training subjects:", length(unique(spirometry$id[trn_idx])), "\n")
#> Training subjects: 407
cat("Test subjects:    ", length(unique(spirometry$id[!trn_idx])), "\n")
#> Test subjects:     102
```

The code below shows the full production fit (M = 200, cv.flag = TRUE).
To keep vignette build time short, pre-computed results are loaded from
the package’s `inst/extdata/` directory when available.

``` r

## Full production fit — stored in inst/extdata/fit_spi.rds
fit_spi <- boostmtree(
  x       = spirometry$features[trn_idx, ],
  tm      = spirometry$time[trn_idx],
  id      = spirometry$id[trn_idx],
  y       = spirometry$y[trn_idx],
  family  = "Continuous",
  M       = 200,
  nu      = 0.05,
  K       = 5,
  cv.flag = TRUE,
  verbose = FALSE
)
```

``` r

print(fit_spi)
#> boostmtree summary
#> model                       : mtree-Pspline learner 
#> fitting mode                : grow 
#> Family                      : Continuous 
#> number of K-terminal nodes  : 5 
#> regularization parameter    : 0.05 
#> sample size                 : 407 
#> number of variables         : 23 
#> number of unique time points: 2510 
#> avg. number of time points  : 18.41 
#> B-spline dimension          : 14 
#> penalization order          : 3 
#> boosting iterations         : 50 
#> optimized number iterations : 50 
#> optimized rho               : 0.7663 
#> optimized phi               : 425.3598 
#> OOB cv RMSE                 : 0.9405
```

### Diagnostics

``` r

plot(fit_spi)
```

The top-left panel shows OOB error (RMSE) vs. boosting iteration. The
dashed vertical line marks the selected `Mopt`. Use this to judge
whether more iterations (`M`) are needed or if the current budget is
sufficient.

### Test-set prediction

``` r

pred_spi <- predict(
  fit_spi,
  x  = spirometry$features[!trn_idx, ],
  tm = spirometry$time[!trn_idx],
  id = spirometry$id[!trn_idx],
  y  = spirometry$y[!trn_idx]
)

print(pred_spi)
#> boostmtree summary
#> model                       : mtree-Pspline learner 
#> fitting mode                : predict 
#> Family                      : Continuous 
#> sample size                 : 102 
#> number of variables         : 23 
#> number of unique time points: 1091 
#> avg. number of time points  : 19.4 
#> optimized number iterations : 19 
#> optimized rho               : 0.7684 
#> optimized phi               : 457.4402 
#> test set RMSE               : 0.8843
```

### Variable importance

``` r

vimp_spi <- vimp.boostmtree(fit_spi)
vimpPlot(vimp_spi)
```

The bar chart separates *main effects* (bars above the axis, labelled in
the upper section) from *time-interaction effects* (bars below the
axis). Variables with both a main effect and a time-interaction effect
are the most interesting: their impact on FEV1% changes over
post-transplant time.

### Marginal dependence

``` r

## Inspect the top-4 most important variables
top4 <- rownames(vimp_spi[[1]])[order(abs(vimp_spi[[1]][, 1]),
                                       decreasing = TRUE)][1:4]
marginalPlot(fit_spi, xvar.names = top4, plot.it = TRUE)
```

------------------------------------------------------------------------

## Atrial Fibrillation: Binary longitudinal response

The `AF` dataset tracks weekly AF presence/absence (binary) for 228
patients enrolled in a surgical ablation trial, with up to 12 months
follow-up.

``` r

data(AF, package = "boostmtree")

cat("Patients:    ", length(unique(AF$id)), "\n")
#> Patients:     228
cat("Total obs:   ", length(AF$y), "\n")
#> Total obs:    7949
cat("Covariates:  ", ncol(AF$feature), "\n")
#> Covariates:   34
cat("AF prevalence: ", round(mean(AF$y == 1) * 100, 1), "%\n")
#> AF prevalence:  52.2 %
```

### Model fit

For binary responses, set `family = "Binary"`. The model fits the
log-odds on the linear predictor scale; predictions are returned as
probabilities.

``` r

set.seed(2024)

subjects_af <- unique(AF$id)
trn_af_subj <- sample(subjects_af, floor(0.75 * length(subjects_af)))
trn_af      <- AF$id %in% trn_af_subj

cat("Training subjects:", length(unique(AF$id[trn_af])), "\n")
#> Training subjects: 171
cat("Test subjects:    ", length(unique(AF$id[!trn_af])), "\n")
#> Test subjects:     57
```

The code below shows the full production fit. Pre-computed results are
loaded from `inst/extdata/fit_af.rds` when available.

``` r

## Full production fit — stored in inst/extdata/fit_af.rds
fit_af <- boostmtree(
  x       = AF$feature[trn_af, ],
  tm      = AF$time[trn_af],
  id      = AF$id[trn_af],
  y       = AF$y[trn_af],
  family  = "Binary",
  M       = 200,
  nu      = 0.05,
  cv.flag = TRUE,
  verbose = FALSE
)
```

``` r

print(fit_af)
#> boostmtree summary
#> model                       : mtree-Pspline learner 
#> fitting mode                : grow 
#> Family                      : Binary 
#> number of K-terminal nodes  : 5 
#> regularization parameter    : 0.05 
#> sample size                 : 171 
#> number of variables         : 34 
#> number of unique time points: 5826 
#> avg. number of time points  : 34.25 
#> B-spline dimension          : 14 
#> penalization order          : 3 
#> boosting iterations         : 50 
#> optimized number iterations : 49 
#> optimized rho               : 0.6377 
#> optimized phi               : 0.2011 
#> OOB cv RMSE                 : 0.4645
```

``` r

plot(fit_af)
```

### Test-set prediction

``` r

pred_af <- predict(
  fit_af,
  x  = AF$feature[!trn_af, ],
  tm = AF$time[!trn_af],
  id = AF$id[!trn_af],
  y  = AF$y[!trn_af]
)

print(pred_af)
#> boostmtree summary
#> model                       : mtree-Pspline learner 
#> fitting mode                : predict 
#> Family                      : Binary 
#> sample size                 : 57 
#> number of variables         : 34 
#> number of unique time points: 2089 
#> avg. number of time points  : 36.72 
#> optimized number iterations : 20 
#> optimized rho               : 0.6572 
#> optimized phi               : 0.2174 
#> test set RMSE               : 0.4683
```

### Variable importance

``` r

vimp_af <- vimp.boostmtree(fit_af)
vimpPlot(vimp_af)
```

### Marginal covariate effects

``` r

marginalPlot(fit_af, xvar.names = colnames(AF$feature)[1:4], plot.it = TRUE)
```

------------------------------------------------------------------------

## Simulation study: comparing models

[`simLong()`](https://ehrlinger.github.io/boostmtree/reference/simLong.md)
supports systematic benchmarking across simulation models. Below we
compare OOB RMSE across the four simulation scenarios for a continuous
response.

``` r

results <- lapply(0:3, function(m) {
  set.seed(m + 100)
  dta <- simLong(n = 50, ntest = 0, N = 4, model = m,
                 family = "Continuous", q = 2)
  fit <- boostmtree(
    x       = dta$dtaL$features,
    tm      = dta$dtaL$time,
    id      = dta$dtaL$id,
    y       = dta$dtaL$y,
    family  = "Continuous",
    M       = 50,
    cv.flag = TRUE,
    verbose = FALSE
  )
  data.frame(
    model    = paste0("Model ", m),
    formula  = dta$f.true,
    Mopt     = fit$Mopt[1],
    OOB_RMSE = round(min(fit$rmse, na.rm = TRUE), 3)
  )
})

do.call(rbind, results)
#>     model                                            formula Mopt OOB_RMSE
#> 1 Model 0                              y ~ g( x1 + x3 + x4 )   50    0.918
#> 2 Model 1               y ~ g( x1 + x3 + x4 + I(time * x2) )   50    0.680
#> 3 Model 2           y ~ g( x1 + x3 + x4 + I(time^2 * x2^2) )   50    0.611
#> 4 Model 3 y ~ g( x1 + x3 + exp(x4) + I(time^2 * x2^2 * x3) )   50    0.597
```

As expected, error increases with model complexity (richer
covariate–time interactions), but boosting adapts `Mopt` accordingly.

------------------------------------------------------------------------

## Tuning guide

| Parameter | Effect | Recommendation |
|:---|:---|:---|
| `M` | Maximum iterations | Start at 200–500; increase if OOB error is still decreasing at `M` |
| `nu` | Shrinkage (learning rate) | Smaller values (0.01–0.05) give better results but need larger `M` |
| `K` | Tree terminal nodes | 3–10; smaller values give smoother fits |
| `nknots` | B-spline knots for time interaction | 5–15; more knots = more flexible time curves |
| `cv.flag` | OOB-based optimal stopping | Always use `TRUE`; safeguards against overfitting |
| `ntree` | Trees per base learner (passed to `randomForestSRC`) | 50–200 |

### Learning rate vs. iterations trade-off

A smaller `nu` makes each update more conservative, reducing overfitting
risk but requiring more iterations to converge. A common strategy is:

``` r

# Initial exploration
fit_fast <- boostmtree(x, tm, id, y, M = 100, nu = 0.1, cv.flag = TRUE)

# Refine with slower learning rate
fit_full <- boostmtree(x, tm, id, y,
                       M   = fit_fast$Mopt * 5,
                       nu  = 0.02,
                       cv.flag = TRUE)
```

------------------------------------------------------------------------

## References

Friedman J.H. (2001). Greedy function approximation: a gradient boosting
machine. *Annals of Statistics*, 29(5): 1189–1232.

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
data. *Machine Learning*, 106(2): 277–305. doi:
[10.1007/s10994-016-5597-1](https://doi.org/10.1007/s10994-016-5597-1)

Gillinov A.M. et al. (2015). Surgical ablation of atrial fibrillation
during mitral valve surgery. *New England Journal of Medicine*, 372(15):
1399–1408.

Mason D.P. et al. (2012). Effect of changes in postoperative spirometry
on survival after lung transplantation. *Journal of Thoracic and
Cardiovascular Surgery*, 144: 197–203.
