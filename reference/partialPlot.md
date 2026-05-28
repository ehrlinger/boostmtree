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
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpkT6Fbi/PartialPlot.pdf
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x2",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpkT6Fbi/PartialPlot.pdf

#partial plot using "x2" as the conditional variable
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",
                      conditional.xvar = "x2", conditional.values = 1,plot.it = TRUE)
#> Plot saved to: /tmp/RtmpkT6Fbi/PartialPlot.pdf
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",
                      conditional.xvar = "x2", conditional.values = 2,plot.it = TRUE)
#> Plot saved to: /tmp/RtmpkT6Fbi/PartialPlot.pdf

##------------------------------------------------------------
## Synthetic example (Response is binary)
## high correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 50, N = 5, rho =.80, model = 2,family = "Binary")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y,family = "Binary",M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========================================                            |  60%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x1",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpkT6Fbi/PartialPlot.pdf
pp.obj <- partialPlot(object = boost.grow, xvar.names = "x2",plot.it = TRUE)
#> Plot saved to: /tmp/RtmpkT6Fbi/PartialPlot.pdf
# }

if (FALSE) { # \dontrun{
##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------
data(spirometry, package = "boostmtree")

#boosting call: cubic B-splines with 15 knots
spr.obj <- boostmtree(spirometry$features, spirometry$time, spirometry$id, spirometry$y,
            family = "Continuous",M = 300, nu = .025, nknots = 15)

#partial plot of double-lung group at 5 years
dltx <- partialPlot(object = spr.obj, xvar.names = "AGE",
                    tm.unq = 5, subset=spr.obj$x$DOUBLE==1,plot.it = TRUE)

#partial plot of single-lung group at 5 years
sltx <- partialPlot(object = spr.obj, xvar.names = "AGE",
                    tm.unq = 5, subset=spr.obj$x$DOUBLE==0,plot.it = TRUE)

#combine the two plots: we use lowess smoothed values
dltx <- dltx$l.obj[[1]]
sltx <- sltx$l.obj[[1]]
plot(range(c(dltx[, 1], sltx[, 1])), range(c(dltx[, -1], sltx[, -1])),
     xlab = "age", ylab = "predicted y (adjusted)", type = "n")
lines(dltx[, 1], dltx[, -1], lty = 1, lwd = 2, col = "red")
lines(sltx[, 1], sltx[, -1], lty = 1, lwd = 2, col = "blue")
legend("topright", legend = c("DLTx", "SLTx"), lty = 1, fill = c(2,4))
} # }
```
