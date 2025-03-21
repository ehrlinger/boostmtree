boostmtree - Boosted multivariate trees for longitudinal data.
===============================================================

<!-- badges: start -->

[![cranlogs](http://cranlogs.r-pkg.org/badges/boostmtree)](http://cranlogs.r-pkg.org/badges/boostmtree)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/boostmtree)](https://cran.r-project.org/package=boostmtree)

  [![R-CMD-check](https://github.com/ehrlinger/boostmtree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ehrlinger/boostmtree/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ehrlinger/boostmtree/graph/badge.svg)](https://app.codecov.io/gh/ehrlinger/boostmtree)
<!-- badges: end -->
  
Multivariate extension of Friedman's (2001) gradient descent boosting method for modeling longitudinal response using multivariate tree base learners.
  Longitudinal response could be continuous, binary, nominal or ordinal.
  Covariate-time interactions are modeled using penalized B-splines
  (P-splines) with estimated adaptive smoothing parameter.

Package Overview
=====
  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.

 1. `boostmtree` - This is the main entry point to the package.  It grows a multivariate tree using user supplied training data.  Trees are grown using the `randomForestSRC` R-package.

 2. `predict.boostmtree` - Used for prediction.  Predicted values are obtained by dropping the user supplied test data down the grow forest. The resulting object has class (`rfsrc`, `predict`).

## Authors
Hemant Ishwaran, Amol Pande and Udaya B. Kogalur

# References

  Friedman J.H. (2001). Greedy function approximation: a gradient boosting machine, `Ann. of Statist.`, 5:1189-1232.
  
  Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U., Blackstone E.H., Ishwaran H. Boosted Multivariate Trees for Longitudinal Data `Mach Learn`,
. 2017 Feb;106(2):277-305.  doi: 10.1007/s10994-016-5597-1.
