# Plot Proximity Profile

Plot the mean longitudinal profile for a randomly selected subject and
its proximity-matched neighbours.

## Usage

``` r
# S3 method for class 'profile.prx'
plot(x, col = NULL, rnd.case = NULL, cut = 0.95, restrictX = TRUE, ...)
```

## Arguments

- x:

  An object of class `profile.prx` returned by
  [`predict.boostmtree`](https://ehrlinger.github.io/boostmtree/reference/predict.boostmtree.md)
  when `proximity = TRUE`.

- col:

  Optional colour vector for matched subjects.

- rnd.case:

  Index of the focal subject. If `NULL` (default), one subject is chosen
  at random.

- cut:

  Quantile threshold (default `0.95`) used to select high-proximity
  neighbours.

- restrictX:

  Logical; if `TRUE` (default) the x-axis is restricted to the time
  range of the focal subject.

- ...:

  Further arguments passed to or from other methods.

## Value

No return value, called for side effects.
