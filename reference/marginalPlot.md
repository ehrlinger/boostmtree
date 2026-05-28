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
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
marginalPlot(boost.grow, "x1",plot.it = TRUE)
#> Plot saved to: /tmp/RtmplzDByJ/MarginalPlot.pdf
marginalPlot(boost.grow, "x2",plot.it = TRUE)
#> Plot saved to: /tmp/RtmplzDByJ/MarginalPlot.pdf

#Plot of all covariates. The plot will be stored as the "MarginalPlot.pdf"
# in the current working directory.
marginalPlot(boost.grow,plot.it = TRUE)
#> Plot saved to: /tmp/RtmplzDByJ/MarginalPlot.pdf


##------------------------------------------------------------
## Synthetic example (Response is binary)
## High correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 50, N = 5, rho =.80, model = 2,family = "Binary")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y, family = "Binary", M = 20)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |=================================================                     |  70%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |====================================================                  |  75%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |========================================================              |  80%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |============================================================          |  85%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |===============================================================       |  90%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================================================================    |  95%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
marginalPlot(boost.grow, "x1",plot.it = TRUE)
#> Plot saved to: /tmp/RtmplzDByJ/MarginalPlot.pdf
marginalPlot(boost.grow, "x2",plot.it = TRUE)
#> Plot saved to: /tmp/RtmplzDByJ/MarginalPlot.pdf

#Plot of all covariates. The plot will be stored as the "MarginalPlot.pdf"
# in the current working directory.
marginalPlot(boost.grow,plot.it = TRUE)
#> Plot saved to: /tmp/RtmplzDByJ/MarginalPlot.pdf
# }

if (FALSE) { # \dontrun{
##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------
data(spirometry, package = "boostmtree")

#boosting call: cubic B-splines with 15 knots
spr.obj <- boostmtree(spirometry$features, spirometry$time, spirometry$id, spirometry$y,
            family = "Continuous",M = 300, nu = .025, nknots = 15)

#marginal plot of double-lung group at 5 years
dltx <- marginalPlot(spr.obj, "AGE", tm.unq = 5, subset = spr.obj$x$DOUBLE==1,plot.it = TRUE)

#marginal plot of single-lung group at 5 years
sltx <- marginalPlot(spr.obj, "AGE", tm.unq = 5, subset = spr.obj$x$DOUBLE==0,plot.it = TRUE)

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
} # }
```
