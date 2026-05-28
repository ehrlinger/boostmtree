# Simulate longitudinal data

Simulates longitudinal data with continuous or binary response from
models with increasing complexity of covariate-time interactions.

## Usage

``` r
simLong(
  n = 100,
  ntest = 0,
  N = 5,
  rho = 0.8,
  type = c("corCompSym", "corAR1", "corSymm", "iid"),
  model = c(0, 1, 2, 3),
  family = c("Continuous", "Binary"),
  phi = 1,
  q = 0,
  ...
)
```

## Arguments

- n:

  Requested training sample size.

- ntest:

  Requested test sample size.

- N:

  Parameter controlling number of time points per subject.

- rho:

  Correlation parameter.

- type:

  Type of correlation matrix.

- model:

  Requested simulation model.

- family:

  Family of response `y`. Use `"Continuous"` or `"Binary"` based on the
  scale of `y`.

- phi:

  Variance of measurement error.

- q:

  Number of zero-signal variables (i.e., variables unrelated to y).

- ...:

  Further arguments passed to or from other methods.

## Value

An invisible list with the following components:

- dtaL:

  List containing the simulated data in the following order: `features`,
  `time`, `id` and `y`.

- dta:

  Simulated data given as a data frame.

- trn:

  Index of `id` values identifying the training data.

- f.true:

  Formula of the simulation model.

## Details

Simulates longitudinal data with 3 main effects and (possibly) a
covariate-time interaction. Complexity of the model is specified using
the option `model`:

1.  *`model=0`:* Linear with no covariate-time interactions.

2.  *`model=1`:* Linear covariate-time interaction.

3.  *`model=2`:* Quadratic time-quadratic covariate interaction.

4.  *`model=3`:* Quadratic time-quadratic two-way covariate interaction.

For details see Pande et al. (2017).

## References

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
data, *Machine Learning*, 106(2): 277–305.

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur

## Examples

``` r

# \donttest{
##------------------------------------------------------------
##  Response is continuous
##----------------------------------------------------------------------------

## set the number of boosting iterations
M <- 20

## simulation 0: only main effects (x1, x3, x4)
dta <- simLong(n = 20, ntest = 5, model = 0, family = "Continuous", q = 5)
trn <- dta$trn
dtaL <- dta$dtaL
dta <- dta$dta
obj.0 <-  boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
          family = "Continuous", M = M)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=====================                                                 |  30%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================                                |  55%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========================================                            |  60%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |====================================================                  |  75%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================================================              |  80%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
pred.0 <- predict(obj.0, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])



##------------------------------------------------------------
##  Response is binary
##----------------------------------------------------------------------------

## set the number of boosting iterations
M <- 20

## simulation 0: only main effects (x1, x3, x4)
dta <- simLong(n = 20, ntest = 5, model = 0, family = "Binary", q = 5)
trn <- dta$trn
dtaL <- dta$dtaL
dta <- dta$dta
obj.0 <-  boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
          family = "Binary", M = M)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============                                                        |  20%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=====================                                                 |  30%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================                                              |  35%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================                                          |  40%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |================================                                      |  45%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#>   |                                                                              |======================================                                |  55%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
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
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%
#> qr.solve failed (VMat): singular matrix 'a' in solve
#> qr.solve failed (VMat): singular matrix 'a' in solve
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
#>   |                                                                              |====================================================                  |  75%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================================================              |  80%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================================================          |  85%
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
#>   |                                                                              |======================================================================| 100%
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
pred.0 <- predict(obj.0, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])
# }
```
