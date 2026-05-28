# Plot Proximity Profile for a `profile.prx` Object

S3 method for plotting proximity-weighted mean profiles from a
`boostmtree` prediction object computed with `proximity = TRUE`. A
randomly selected (or user-specified) subject is used as the reference
case. Subjects whose proximity score to the reference case exceeds the
`cut` quantile are treated as neighbours; their smoothed fitted
trajectories (dashed lines) and a proximity-weighted average trajectory
(solid black line) are overlaid on the reference subject's smoothed
trajectory (solid blue line).

## Usage

``` r
# S3 method for class 'profile.prx'
plot(x, col = NULL, rnd.case = NULL, cut = 0.95, restrictX = TRUE, ...)
```

## Arguments

- x:

  An object of class `profile.prx` — a `boostmtree` predict object that
  contains a `proximity` matrix (i.e., produced with
  `proximity = TRUE`).

- col:

  Optional colour vector. Must be indexed by *original subject index* —
  colours are looked up as `col[i]` where `i` is the subject's row
  number in the proximity matrix. A safe choice is a vector of length
  `nrow(x$proximity)`. Defaults to `1` (black) for every matched
  neighbour.

- rnd.case:

  Integer index of the reference subject. If `NULL` (default) a subject
  is chosen at random.

- cut:

  Quantile threshold (0–1) for defining neighbours based on proximity
  score. Defaults to `0.95`.

- restrictX:

  Logical. If `TRUE` (default) the x-axis is restricted to the
  observation times of the reference subject.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (currently
  unused).

## Value

Invisibly returns the covariate matrix rows corresponding to the matched
neighbours.

## See also

[`predict.boostmtree`](https://ehrlinger.github.io/boostmtree/reference/predict.boostmtree.md)
