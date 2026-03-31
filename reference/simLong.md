# Simulate longitudinal data

Simulates longitudinal data with continuous or binary response from
models with increasing complexity of covariate-time interactions.

## Usage

``` r
simLong(n,
        ntest = 0,
        N = 5,
        rho = 0.8,
        type = c("corCompSym", "corAR1", "corSymm", "iid"),
        model = c(0, 1, 2, 3),
        family = c("Continuous","Binary"),
        phi = 1,
        q = 0,
        ...)
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

  Family of response `y`. Use any one from {"Continuous", "Binary"}
  based on the scale of `y`.

- phi:

  Variance of measurement error.

- q:

  Number of zero-signal variables (i.e., variables unrelated to y).

- ...:

  Further arguments passed to or from other methods.

## Details

Simulates longitudinal data with 3 main effects and (possibly) a
covariate-time interaction. Complexity of the model is specified using
the option `model`:

1.  *`model=0`:* Linear with no covariate-time interactions.

2.  *`model=1`:* Linear covariate-time interaction.

3.  *`model=2`:* Quadratic time-quadratic covariate interaction.

4.  *`model=3`:* Quadratic time-quadratic two-way covariate interaction.

For details see Pande et al. (2017).

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

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur

## References

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
data, *Machine Learning*, 106(2): 277–305.

## Examples

``` r
if (FALSE) { # \dontrun{
##------------------------------------------------------------
##  Response is continuous
##----------------------------------------------------------------------------

## set the number of boosting iterations
M <- 500

## simulation 0: only main effects (x1, x3, x4)
dta <- simLong(n = 100, ntest = 100, model = 0, family = "Continuous", q = 5)
trn <- dta$trn
dtaL <- dta$dtaL
dta <- dta$dta
obj.0 <-  boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn], 
          family = "Continuous", M = M)
pred.0 <- predict(obj.0, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])



##------------------------------------------------------------
##  Response is binary
##----------------------------------------------------------------------------

## set the number of boosting iterations
M <- 500

## simulation 0: only main effects (x1, x3, x4)
dta <- simLong(n = 100, ntest = 100, model = 0, family = "Binary", q = 5)
trn <- dta$trn
dtaL <- dta$dtaL
dta <- dta$dta
obj.0 <-  boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn], 
          family = "Binary", M = M)
pred.0 <- predict(obj.0, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])
} # }
```
