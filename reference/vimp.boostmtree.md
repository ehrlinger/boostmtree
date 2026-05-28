# Variable Importance

Calculate VIMP score for each of the individual covariates or a joint
VIMP of multiple covariates.

## Usage

``` r
vimp.boostmtree(object, x.names = NULL, joint = FALSE)
```

## Arguments

- object:

  A boosting object of class `(boostmtree, grow)` or class
  `(boostmtree, predict)`.

- x.names:

  Names of the x-variables for which VIMP is requested. If `NULL`, VIMP
  is calculated for all covariates.

- joint:

  Logical. If `FALSE` (default), individual VIMP is returned for each
  covariate in `x.names`. If `TRUE`, a single joint VIMP is computed for
  all covariates combined.

## Value

A list with three components:

- vimp.main:

  Matrix of main-effect VIMP scores (rows = variables, columns =
  response classes).

- vimp.int:

  Matrix of covariate-time interaction VIMP scores.

- vimp.time:

  Numeric vector of pure time-effect VIMP scores.

For cross-sectional (univariate) data, only `vimp.main` is populated.

## Details

Variable Importance (VIMP) is calculated for each of the covariates
individually or a joint VIMP is calculated for all the covariates
specified in `x.names`.

## References

Friedman J.H. Greedy function approximation: a gradient boosting
machine, *Ann. of Statist.*, 5:1189-1232, 2001.

## See also

[`vimpPlot`](https://ehrlinger.github.io/boostmtree/reference/vimpPlot.md),
[`boostmtree`](https://ehrlinger.github.io/boostmtree/reference/boostmtree.md),
[`predict.boostmtree`](https://ehrlinger.github.io/boostmtree/reference/predict.boostmtree.md)

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur

## Examples

``` r

# \donttest{
##------------------------------------------------------------
## Synthetic example (Response is continuous)
## VIMP is based on in-sample CV using out of bag data
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 20, N = 5, rho =.80, model = 2,family = "Continuous")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y,
              family = "Continuous", M = 20,cv.flag = TRUE)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================                                          |  40%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================                                |  55%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========================================                            |  60%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |====================================================                  |  75%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%
vimp.grow <- vimp.boostmtree(object = boost.grow,x.names=c("x1","x2"),joint = FALSE)
vimp.joint.grow <- vimp.boostmtree(object = boost.grow,x.names=c("x1","x2"),joint = TRUE)

##------------------------------------------------------------
## Synthetic example (Response is continuous)
## VIMP is based on test data
##-------------------------------------------------------------
#simulate the data
dtaO <- simLong(n = 20, ntest = 10, N = 5, rho =.80, model = 2, family = "Continuous")

## save the data as both a list and data frame
dtaL <- dtaO$dtaL
dta <- dtaO$dta

## get the training data
trn <- dtaO$trn

#basic boosting call
boost.grow <- boostmtree(dtaL$features[trn,], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
              family = "Continuous", M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========================================                            |  60%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%
boost.pred <- predict(boost.grow,dtaL$features[-trn,], dtaL$time[-trn], dtaL$id[-trn],
              dtaL$y[-trn])
vimp.pred <- vimp.boostmtree(object = boost.pred,x.names=c("x1","x2"),joint = FALSE)
vimp.joint.pred <- vimp.boostmtree(object = boost.pred,x.names=c("x1","x2"),joint = TRUE)

# }
```
