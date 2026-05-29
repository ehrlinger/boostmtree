# Marginal plot analysis

Marginal plot of x against the unadjusted predicted y. This is mainly
used to obtain marginal relationships between x and the unadjusted
predicted y. Marginal plots have a faster execution compared to partial
plots (Friedman, 2001).

## Usage

``` r
marginalPlot(
  object,
  xvar.names,
  tm.unq,
  subset,
  plot.it = FALSE,
  path_saveplot = NULL,
  Verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A boosting object of class `(boostmtree, grow)`.

- xvar.names:

  Names of the x-variables to be used. By default, all variables are
  plotted.

- tm.unq:

  Unique time points used for the plots of x against y. By default, the
  deciles of the observed time values are used.

- subset:

  Vector indicating which rows of the x-data to be used for the
  analysis. The default is to use the entire data.

- plot.it:

  Should plots be displayed? If `xvar.names` is a vector with more than
  one variable name, then instead of displaying, plot is stored as
  "MarginalPlot.pdf" in the location specified by `path_saveplot`.

- path_saveplot:

  Provide the location where plot should be saved. By default the plot
  will be saved at temporary folder.

- Verbose:

  Display the path where the plot is saved?

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly returns a list with components `p.obj` (marginal effect
estimates), `l.obj` (lowess smoothed longitudinal estimates, or `NULL`
if `plot.it = FALSE`), and `time` (time points used for evaluation).

## Details

Marginal plot of x values specified by `xvar.names` against the
unadjusted predicted y-values over a set of time points specified by
`tm.unq`. Analysis can be restricted to a subset of the data using
`subset`.

## References

Friedman J.H. Greedy function approximation: a gradient boosting
machine, *Ann. of Statist.*, 5:1189-1232, 2001.

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur

## Examples

``` r

# \donttest{
##------------------------------------------------------------
## Synthetic example (Response is continuous)
## High correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 50, N = 5, rho =.80, model = 2,family = "Continuous")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y, family = "Continuous", M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
marginalPlot(boost.grow, "x1",plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf
marginalPlot(boost.grow, "x2",plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf

#Plot of all covariates. The plot will be stored as the "MarginalPlot.pdf"
# in the current working directory.
marginalPlot(boost.grow,plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf


##------------------------------------------------------------
## Synthetic example (Response is binary)
## High correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 50, N = 5, rho =.80, model = 2,family = "Binary")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y, family = "Binary", M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================                                              |  35%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================================================              |  80%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
marginalPlot(boost.grow, "x1",plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf
marginalPlot(boost.grow, "x2",plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf

#Plot of all covariates. The plot will be stored as the "MarginalPlot.pdf"
# in the current working directory.
marginalPlot(boost.grow,plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf
# }

# \donttest{
##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------
data(spirometry, package = "boostmtree")

#boosting call: cubic B-splines with 15 knots
spr.obj <- boostmtree(spirometry$features, spirometry$time, spirometry$id, spirometry$y,
            family = "Continuous",M = 300, nu = .025, nknots = 15)
#>   |                                                                              |                                                                      |   0%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |                                                                      |   1%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=                                                                     |   1%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=                                                                     |   2%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==                                                                    |   2%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==                                                                    |   3%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===                                                                   |   4%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===                                                                   |   5%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====                                                                  |   5%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====                                                                  |   6%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====                                                                 |   7%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====                                                                 |   8%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======                                                                |   8%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======                                                                |   9%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======                                                               |   9%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======                                                               |  10%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======                                                               |  11%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========                                                              |  11%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========                                                              |  12%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========                                                             |  12%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========                                                             |  13%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========                                                            |  14%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========                                                            |  15%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========                                                           |  15%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========                                                           |  16%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============                                                          |  17%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============                                                          |  18%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============                                                         |  18%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============                                                         |  19%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============                                                        |  19%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============                                                        |  20%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============                                                        |  21%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============                                                       |  21%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============                                                       |  22%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================                                                      |  22%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================                                                      |  23%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================                                                     |  24%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================                                                     |  25%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================                                                    |  25%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================                                                    |  26%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================                                                   |  27%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================                                                   |  28%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================                                                  |  28%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================                                                  |  29%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================                                                 |  29%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================                                                 |  30%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================                                                 |  31%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================                                                |  31%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================                                                |  32%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======================                                               |  32%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======================                                               |  33%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================                                              |  34%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================                                              |  35%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================                                             |  35%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================                                             |  36%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========================                                            |  37%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========================                                            |  38%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========================                                           |  38%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========================                                           |  39%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================                                          |  39%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================                                          |  40%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================                                          |  41%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============================                                         |  41%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============================                                         |  42%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============================                                        |  42%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============================                                        |  43%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================                                       |  44%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================                                       |  45%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================================                                      |  45%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================================                                      |  46%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================================                                     |  47%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================================                                     |  48%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================================                                    |  48%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================================                                    |  49%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================                                   |  49%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================                                   |  50%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================                                   |  51%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================================                                  |  51%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================================                                  |  52%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================================                                 |  52%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================================                                 |  53%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================================                                |  54%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================================                                |  55%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======================================                               |  55%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======================================                               |  56%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================================                              |  57%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================================                              |  58%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================================                             |  58%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================================                             |  59%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========================================                            |  59%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========================================                            |  60%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========================================                            |  61%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========================================                           |  61%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========================================                           |  62%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================================                          |  62%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================================                          |  63%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============================================                         |  64%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============================================                         |  65%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============================================                        |  65%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============================================                        |  66%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================================                       |  67%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================================                       |  68%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================================================                      |  68%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================================================                      |  69%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================================================                     |  69%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================================================                     |  70%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================================================                     |  71%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================================================                    |  71%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================================================                    |  72%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================================                   |  72%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================================                   |  73%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================================================                  |  74%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================================================                  |  75%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================================================                 |  75%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================================================                 |  76%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================================================                |  77%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================================================                |  78%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======================================================               |  78%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======================================================               |  79%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================================================              |  79%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================================================              |  80%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================================================              |  81%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================================================             |  81%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================================================             |  82%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========================================================            |  82%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========================================================            |  83%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========================================================           |  84%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========================================================           |  85%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================================================          |  86%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============================================================         |  87%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============================================================         |  88%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============================================================        |  88%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============================================================        |  89%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================================================       |  89%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================================================       |  90%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================================================       |  91%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================================================================      |  91%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================================================================      |  92%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================================================================     |  92%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================================================================     |  93%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================================================================    |  94%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================================================================    |  95%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================================================   |  95%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================================================   |  96%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================================================================  |  97%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================================================================  |  98%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================================================== |  98%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================================================================== |  99%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================================================================|  99%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================================================================| 100%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> 
#> gls (full model) failed: computed "gls" fit is singular, rank 25

#marginal plot of double-lung group at 5 years
dltx <- marginalPlot(spr.obj, "AGE", tm.unq = 5, subset = spr.obj$x$DOUBLE==1,plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf

#marginal plot of single-lung group at 5 years
sltx <- marginalPlot(spr.obj, "AGE", tm.unq = 5, subset = spr.obj$x$DOUBLE==0,plot.it = TRUE)
#> Plot saved to: /tmp/Rtmp6dBw9y/MarginalPlot.pdf

#combine the two plots
dltx <- dltx[[2]][[1]]
sltx <- sltx[[2]][[1]]
plot(range(c(dltx[[1]][, 1], sltx[[1]][, 1])), range(c(dltx[[1]][, -1], sltx[[1]][, -1])),
     xlab = "age", ylab = "predicted y", type = "n")
lines(dltx[[1]][, 1][order(dltx[[1]][, 1]) ], dltx[[1]][, -1][order(dltx[[1]][, 1]) ],
      lty = 1, lwd = 2, col = "red")
lines(sltx[[1]][, 1][order(sltx[[1]][, 1]) ], sltx[[1]][, -1][order(sltx[[1]][, 1]) ],
      lty = 1, lwd = 2, col = "blue")
legend("topright", legend = c("DLTx", "SLTx"), lty = 1, fill = c(2,4))

# }
```
