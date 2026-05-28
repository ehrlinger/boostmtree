# Prediction for Boosted multivariate trees for longitudinal data.

Obtain predicted values. Also returns test-set performance if the test
data contains y-outcomes.

## Usage

``` r
# S3 method for class 'boostmtree'
predict(object, x, tm, id, y, M, eps = 1e-05, useCVflag = FALSE, ...)
```

## Arguments

- object:

  A boosting object of class `(boostmtree, grow)`.

- x:

  Data frame (or matrix) containing test set x-values. Rows must be
  duplicated to match the number of time points for an individual. If
  missing, the training x values are used and `tm`, `id` and `y` are not
  required and no performance values are returned.

- tm:

  Time values for each test set individual with one entry for each row
  of `x`. Optional, but if missing, the set of unique time values from
  the training values are used for each individual and no test-set
  performance values are returned.

- id:

  Unique subject identifier, one entry for each row in `x`. Optional,
  but if missing, each individual is assumed to have a full time-profile
  specified by the unique time values from the training data.

- y:

  Test set y-values, with one entry for each row in `x`.

- M:

  Fixed value for the boosting step number. Leave this empty to
  determine the optimized value obtained by minimizing test-set error.

- eps:

  Tolerance value used for determining the optimal `M`. For experts
  only.

- useCVflag:

  Should the predicted value be based on the estimate derived from oob
  sample?

- ...:

  Further arguments passed to or from other methods.

## Value

An object of class `(boostmtree, predict)`, which is a list with the
following components:

- boost.obj:

  The original boosting object.

- x:

  The test x-values, but with only one row per individual (i.e.
  duplicated rows are removed).

- time:

  List with each component containing the time points for a given test
  individual.

- id:

  Sorted subject identifier.

- y:

  List containing the test y-values.

- Y:

  y-values, in the list-format, where nominal or ordinal Response is
  converted into the binary response.

- family:

  Family of `y`.

- ymean:

  Overall mean of y-values for all individuals. If `family` = "Binary",
  "Nominal" or "Ordinal", `ymean` = 0.

- ysd:

  Overall standard deviation of y-values for all individuals. If
  `family` = "Binary", "Nominal" or "Ordinal", `ysd` = 1.

- xvar.names:

  X-variable names.

- K:

  Number of terminal nodes.

- n:

  Total number of subjects.

- ni:

  Number of repeated measures for each subject.

- n.Q:

  Number of class labels for non-continuous response.

- Q_set:

  Class labels for the non-continuous response.

- y.unq:

  Unique y values for the non-continous response.

- nu:

  Boosting regularization parameter.

- D:

  Design matrix for each subject.

- df.D:

  Number of columns of `D`.

- time.unq:

  Vector of the unique time points.

- baselearner:

  List of length *M* containing the base learners.

- gamma:

  List of length *M*, with each component containing the boosted tree
  fitted values.

- membership:

  List of length *M*, with each component containing the terminal node
  membership for a given boosting iteration.

- mu:

  Estimated mean profile at the optimized `M`.

- Prob_class:

  For family == "Ordinal", this provides individual probabilty rather
  than cumulative probabilty.

- muhat:

  Extrapolated mean profile to all unique time points evaluated at the
  the optimized `M`.

- Prob_hat_class:

  Extrapolated `Prob_class` to all unique time points evaluated at the
  the optimized `M`.

- err.rate:

  Test set standardized l1-error and RMSE.

- rmse:

  Test set standardized RMSE at the optimized `M`.

- Mopt:

  The optimized `M`.

## Details

The predicted time profile and performance values are obtained for test
data from the boosted object grown on the training data.

R-side parallel processing is implemented by replacing the R function
`lapply` with `mclapply` found in the parallel package. You can set the
number of cores accessed by `mclapply` by issuing the command
`options(mc.cores = x)`, where `x` is the number of cores. As an
example, issuing the following options command uses all available cores:

`options(mc.cores=detectCores())`

However, this can create high RAM usage, especially when using function
`partialPlot` which calls the `predict` function.

Note that all performance values (for example prediction error) are
standardized by the overall y-standard deviation. Thus, reported RMSE
(root-mean-squared-error) is actually standardized RMSE. Values are
reported at the optimal stopping time.

## References

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
data, *Machine Learning*, 106(2): 277–305.

## See also

[`plot.boostmtree`](https://ehrlinger.github.io/boostmtree/reference/plot.boostmtree.md),
[`print.boostmtree`](https://ehrlinger.github.io/boostmtree/reference/print.boostmtree.md)

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur

## Examples

``` r

# \donttest{
##------------------------------------------------------------
## Synthetic example (Response is continuous)
##
##  High correlation, quadratic time with quadratic interaction
##  largish number of noisy variables
##
##  Illustrates how modified gradient improves performance
##  also compares performance to ideal and well specified linear models
##----------------------------------------------------------------------------

## simulate the data
## simulation 2: main effects (x1, x3, x4), quad-time-interaction (x2)
dtaO <- simLong(n = 20, ntest = 10, model = 2, family = "Continuous", q = 5)

## save the data as both a list and data frame
dtaL <- dtaO$dtaL
dta <- dtaO$dta

## get the training data
trn <- dtaO$trn

## save formulas for linear model comparisons
f.true <- dtaO$f.true
f.linr <- "y~g( x1+x2+x3+x4+x1*time+x2*time+x3*time+x4*time )"


## modified tree gradient (default)
o.1 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn],dtaL$y[trn],
       family = "Continuous",M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
p.1 <- predict(o.1, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])

## non-modified tree gradient (nmtg)
o.2 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
       family = "Continuous",M = 20, mod.grad = FALSE)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========================================                            |  60%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
p.2 <- predict(o.2, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])

## set rho = 0
o.3 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
       family = "Continuous",M = 20, rho = 0)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================                                          |  40%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%
p.3 <- predict(o.3, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])


##rmse values compared to generalized least squares (GLS)
##for true model and well specified linear models (LM)
cat("true LM           :", boostmtree:::gls.rmse(f.true,dta,trn),"\n")
#> true LM           : NA 
cat("well specified LM :", boostmtree:::gls.rmse(f.linr,dta,trn),"\n")
#> well specified LM : NA 
cat("boostmtree        :", p.1$rmse,"\n")
#> boostmtree        : 0.7273789 
cat("boostmtree  (nmtg):", p.2$rmse,"\n")
#> boostmtree  (nmtg): 0.8536049 
cat("boostmtree (rho=0):", p.3$rmse,"\n")
#> boostmtree (rho=0): 0.8134829 

##predicted value plots
plot(p.1)
#> Plot saved to: /tmp/RtmplzDByJ/boostmtree_plot.pdf
plot(p.2)
#> Plot saved to: /tmp/RtmplzDByJ/boostmtree_plot.pdf
plot(p.3)
#> Plot saved to: /tmp/RtmplzDByJ/boostmtree_plot.pdf



##------------------------------------------------------------
## Synthetic example (Response is binary)
##
##  High correlation, quadratic time with quadratic interaction
##  largish number of noisy variables
##----------------------------------------------------------------------------

## simulate the data
## simulation 2: main effects (x1, x3, x4), quad-time-interaction (x2)
dtaO <- simLong(n = 20, ntest = 10, model = 2, family = "Binary", q = 5)

## save the data as both a list and data frame
dtaL <- dtaO$dtaL
dta <- dtaO$dta

## get the training data
trn <- dtaO$trn

## save formulas for linear model comparisons
f.true <- dtaO$f.true
f.linr <- "y~g( x1+x2+x3+x4+x1*time+x2*time+x3*time+x4*time )"


## modified tree gradient (default)
o.1 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn],dtaL$y[trn],
       family = "Binary",M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
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
#>   |                                                                              |==========                                                            |  15%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============                                                        |  20%
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
#>   |                                                                              |=====================                                                 |  30%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================                                              |  35%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================                                          |  40%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================                                |  55%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========================================                            |  60%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |====================================================                  |  75%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================================================              |  80%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
p.1 <- predict(o.1, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])

# }
```
