# Print Summary Output

Print summary output from the boosting analysis.

## Usage

``` r
# S3 method for class 'boostmtree'
print(x, ...)
```

## Arguments

- x:

  An object of class `(boostmtree, grow)` or `(boostmtree, predict)`.

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly returns the `boostmtree` object `x`, following the convention
for print methods.

## References

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone
E.H., Ishwaran H. (2017). Boosted multivariate trees for longitudinal
data, *Machine Learning*, 106(2): 277–305.

## Author

Hemant Ishwaran, Amol Pande and Udaya B. Kogalur
