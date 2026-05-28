# Plot Summary Analysis

Plot summary analysis of the boosting analysis.

## Usage

``` r
# S3 method for class 'boostmtree'
plot(x, use.rmse = TRUE, path_saveplot = NULL, Verbose = TRUE, ...)
```

## Arguments

- x:

  An object of class `(boostmtree, grow)` or `(boostmtree, predict)`.

- use.rmse:

  Report performance values in terms of standardized
  root-mean-squared-error (RMSE) or mean-squared-error (MSE)? Default is
  standardized RMSE.

- path_saveplot:

  Provide the location where plot should be saved. By default the plot
  will be saved at temporary folder.

- Verbose:

  Display the path where the plot is saved?

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly returns the `boostmtree` object `x`.

## Details

Plot summary output, including predicted values and residuals. Also
plots various parameters against the number of boosting iterations.

## References

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
data, *Machine Learning*, 106(2): 277–305.

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur
