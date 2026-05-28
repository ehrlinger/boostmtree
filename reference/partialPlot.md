# Partial plot analysis

Partial dependence plot of x against adjusted predicted y.

## Usage

``` r
partialPlot(
  object,
  M = NULL,
  xvar.names,
  tm.unq,
  xvar.unq = NULL,
  npts = 25,
  subset,
  prob.class = FALSE,
  conditional.xvars = NULL,
  conditional.values = NULL,
  plot.it = FALSE,
  Variable_Factor = FALSE,
  path_saveplot = NULL,
  Verbose = TRUE,
  useCVflag = FALSE,
  ...
)
```

## Arguments

- object:

  A boosting object of class `(boostmtree, grow)`.

- M:

  Fixed value for the boosting step number. If NULL, then use Mopt if it
  is available from the object, else use M

- xvar.names:

  Names of the x-variables to be used. By default, all variables are
  plotted.

- tm.unq:

  Unique time points used for the plots of x against y. By default, the
  deciles of the observed time values are used.

- xvar.unq:

  Unique values used for the partial plot. Default is NULL in which case
  unique values are obtained uniformaly based on the range of variable.
  Values must be provided using list with same length as lenght of
  `xvar.names`.

- npts:

  Maximum number of points used for x. Reduce this value if plots are
  slow.

- subset:

  Vector indicating which rows of the x-data to be used for the
  analysis. The default is to use the entire data.

- prob.class:

  In case of ordinal response, use class probability rather than
  cumulative probability.

- conditional.xvars:

  Vector of character values indicating names of the x-variables to be
  used for further conditioning (adjusting) the predicted y values.
  Variable names should be different from `xvar.names`.

- conditional.values:

  Vector of values taken by the variables from `conditional.xvars`. The
  length of the vector should be same as the length of the vector for
  `conditional.xvars`, which means only one value per conditional
  variable.

- plot.it:

  Should plots be displayed?

- Variable_Factor:

  Default is FALSE. Use TRUE if the variable specified in `xvar.names`
  is a factor.

- path_saveplot:

  Provide the location where plot should be saved. By default the plot
  will be saved at temporary folder.

- Verbose:

  Display the path where the plot is saved?

- useCVflag:

  Should the predicted value be based on the estimate derived from oob
  sample?

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly returns a list with components `p.obj` (partial effect
estimates), `l.obj` (lowess smoothed partial plots, or `NULL` if
`plot.it = FALSE`), and `time` (time points used for evaluation).

## Details

Partial dependence plot (Friedman, 2001) of x values specified by
`xvar.names` against the adjusted predicted y-values over a set of time
points specified by `tm.unq`. Analysis can be restricted to a subset of
the data using `subset`. Further conditioning can be imposed using
`conditional.xvars`.

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
## high correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 50, N = 5, rho =.80, model = 2,family = "Continuous")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y,family = "Continuous",M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x2",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf

#partial plot using "x2" as the conditional variable
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",
                      conditional.xvar = "x2", conditional.values = 1,plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",
                      conditional.xvar = "x2", conditional.values = 2,plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf

##------------------------------------------------------------
## Synthetic example (Response is binary)
## high correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 50, N = 5, rho =.80, model = 2,family = "Binary")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y,family = "Binary",M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x2",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=                                                                     |   1%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=                                                                     |   2%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==                                                                    |   2%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==                                                                    |   3%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===                                                                   |   4%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===                                                                   |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====                                                                  |   6%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====                                                                 |   7%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====                                                                 |   8%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======                                                                |   8%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======                                                                |   9%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======                                                               |   9%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======                                                               |  10%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=======                                                               |  11%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========                                                              |  11%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========                                                              |  12%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========                                                             |  12%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========                                                             |  13%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==========                                                            |  14%
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============                                                          |  17%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============                                                          |  18%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============                                                         |  18%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============                                                         |  19%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============                                                        |  19%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==============                                                        |  21%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============                                                       |  21%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============                                                       |  22%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=================                                                     |  25%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================                                                    |  25%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |==================                                                    |  26%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================                                                   |  27%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===================                                                   |  28%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================                                                  |  28%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================                                                |  31%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
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
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================                                              |  35%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================                                             |  35%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=========================                                             |  36%
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
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
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |====================================                                  |  52%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================================                                 |  52%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=====================================                                 |  53%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
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
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |========================================                              |  57%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================================                       |  68%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |================================================                      |  68%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===========================================================           |  85%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================================================          |  85%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |============================================================          |  86%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |=============================================================         |  87%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
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
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |===============================================================       |  90%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
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
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#>   |                                                                              |======================================================================| 100%
#> gls (full model) failed: computed "gls" fit is singular, rank 25
#> 
#> gls (full model) failed: computed "gls" fit is singular, rank 25

#partial plot of double-lung group at 5 years
dltx <- partialPlot(object = spr.obj, xvar.names = "AGE",
                    tm.unq = 5, subset=spr.obj$x$DOUBLE==1,plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf

#partial plot of single-lung group at 5 years
sltx <- partialPlot(object = spr.obj, xvar.names = "AGE",
                    tm.unq = 5, subset=spr.obj$x$DOUBLE==0,plot.it = TRUE)
#> Plot saved to: /tmp/RtmpAxyv8H/PartialPlot.pdf

#combine the two plots: we use lowess smoothed values
dltx <- dltx$l.obj[[1]]
sltx <- sltx$l.obj[[1]]
plot(range(c(dltx[, 1], sltx[, 1])), range(c(dltx[, -1], sltx[, -1])),
     xlab = "age", ylab = "predicted y (adjusted)", type = "n")
lines(dltx[, 1], dltx[, -1], lty = 1, lwd = 2, col = "red")
lines(sltx[, 1], sltx[, -1], lty = 1, lwd = 2, col = "blue")
legend("topright", legend = c("DLTx", "SLTx"), lty = 1, fill = c(2,4))

# }
```
