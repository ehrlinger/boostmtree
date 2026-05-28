# Variable Importance (VIMP) plot

Barplot displaying VIMP. If the analysis is for the univariate case,
VIMP is displayed above the x-axis. If the analysis is for the
longitudinal case, VIMP for covariates (main effects) are shown above
the x-axis while VIMP for covariate-time interactions (time interaction
effects) are shown below the x-axis. In either case, negative vimp value
is set to zero.

## Usage

``` r
vimpPlot(
  vimp,
  Q_set = NULL,
  Time_Interaction = TRUE,
  xvar.names = NULL,
  cex.xlab = NULL,
  ymaxlim = 0,
  ymaxtimelim = 0,
  subhead.cexval = 1,
  yaxishead = NULL,
  xaxishead = NULL,
  main = "Variable Importance (%)",
  col = grey(0.8),
  cex.lab = 1.5,
  subhead.labels = c("Time-Interactions Effects", "Main Effects"),
  ylbl = FALSE,
  seplim = NULL,
  eps = 0.1,
  Width_Bar = 1,
  path_saveplot = NULL,
  Verbose = TRUE
)
```

## Arguments

- vimp:

  VIMP values.

- Q_set:

  Provide names for various levels of nominal or ordinal response.

- Time_Interaction:

  Whether VIMP is estimated from a longitudinal data, in which case VIMP
  is available for covariate and covariate-time interaction. Default is
  TRUE. If FALSE, VIMP is assumed to be estimated from a cross-sectional
  data.

- xvar.names:

  Names of the covariates. If NULL, names are assigned as x1, x2,...,xp.

- cex.xlab:

  Magnification of the names of the covariates above (and below) the
  barplot.

- ymaxlim:

  By default, we use the range of the vimp values for the covariates for
  the ylim. If one wants to extend the ylim, add the amount with which
  the ylim will extend above.

- ymaxtimelim:

  By default, we use the range of the vimp values for the
  covariates-time for the ylim. If one wants to extend the ylim, add the
  amount with which the ylim will extend below. Argument only works for
  the longitudinal setting.

- subhead.cexval:

  Magnification of the `subhead.labels`. Argument only works for the
  longitudinal setting.

- yaxishead:

  This represent a vector with two values which are points on the
  y-axis. Corresponding to the values, the lables for `subhead.labels`
  is shown. First argument corresponds to covariate-time interaction,
  whereas second argument is for the main effect. Argument only works
  for the longitudinal setting.

- xaxishead:

  This represent a vector with two values which are points on the
  x-axis. Corresponding to the values, the lables for `subhead.labels`
  is shown. First argument corresponds to covariate-time interaction,
  whereas second argument is for the main effect. Argument only works
  for the longitudinal setting.

- main:

  Main title for the plot.

- col:

  Color of the plot.

- cex.lab:

  Magnification of the x and y lables.

- subhead.labels:

  Labels corresponding to the plot. Default is "Time-Interactions
  Effects" for the barplot below x-axis, and "Main Effects" for the
  barplot above x-axis.

- ylbl:

  Should labels for the sub-headings be shown on left side of the
  y-axis.

- seplim:

  if `ylbl` is `TRUE`, the distance between the lables of the
  sub-headings.

- eps:

  Amount of gap between the top of the barplot and variable names.

- Width_Bar:

  Width of the barplot.

- path_saveplot:

  Provide the location where plot should be saved. By default the plot
  will be saved at temporary folder.

- Verbose:

  Display the path where the plot is saved?

## Value

Invisibly returns the vimp list used for plotting.

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur

## Examples

``` r

# \donttest{
##------------------------------------------------------------
## Synthetic example
## high correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 20, N = 5, rho =.80, model = 2,family = "Continuous")$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y,
              family = "Continuous",M = 20, cv.flag = TRUE)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%
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
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%
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
#>   |                                                                              |==========================================                            |  60%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#>   |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%
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
#>   |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%
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
#>   |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
#> qr.solve failed (Hessian NR): singular matrix 'a' in solve
vimp.grow <- vimp.boostmtree(object = boost.grow)

# VIMP plot
vimpPlot(vimp = vimp.grow, ymaxlim = 20, ymaxtimelim = 20,
         xaxishead = c(3,3), yaxishead = c(65,65),
         cex.xlab = 1, subhead.cexval = 1.2)
#> Warning: length(labels) > max(length(x), length(y));
#> 'labels' truncated to length 4
#> Plot saved to: /tmp/RtmpAxyv8H/VIMPplot.pdf
# }
```
