\name{predict.boostmtree}
\alias{predict.boostmtree}
\title{Prediction for Boosted multivariate trees for longitudinal data.}
\description{
  Obtain predicted values.  Also returns test-set performance
  if the test data contains y-outcomes.
}
\usage{\method{predict}{boostmtree}(object,
 x,
 tm,
 id,
 y,
 M,
 importance = TRUE,
 eps = 1e-5,
 ...)
}
\arguments{
  \item{object}{A boosting object of class \code{(boostmtree, grow)}.}
  \item{x}{Data frame (or matrix) containing test set x-values.  Rows
    must be duplicated to match the number of time points for an
    individual. If missing, the training x values are used and
    \code{tm}, \code{id} and \code{y} are not required and no
    performance values are returned.}
  \item{tm}{Time values for each test set individual with one entry for
    each row of \code{x}.  Optional, but if missing, the set of unique time
    values from the training values are used for each individual and no
    test-set performance values are returned.}
  \item{id}{Unique subject identifier, one entry for each row in
    \code{x}.  Optional, but if missing, each individual is assumed to
    have a full time-profile specified by the unique time values from
    the training data.}
  \item{y}{Test set y-values, with one entry for each row  in \code{x}.}
  \item{M}{Fixed value for the boosting step number.  Leave this empty 
    to determine the optimized value obtained by minimizing test-set error.}
  \item{importance}{Calculate importance values for the x-variables?}
  \item{eps}{Tolerance value used for determining the optimal \code{M}.
    For experts only.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  The predicted time profile and performance values are obtained for
  test data from the boosted object grown on the training data.

  R-side parallel processing is implemented by replacing the R function
  \command{lapply} with \command{mclapply} found in the \pkg{parallel}
  package.  You can set the number of cores accessed by
  \command{mclapply} by issuing the command \command{options(mc.cores =
  x)}, where \command{x} is the number of cores.  As an example, issuing
  the following options command uses all available cores:
  
  \command{options(mc.cores=detectCores())}
  
  However, this can create high RAM usage, especially when using
  function \command{partialPlot} which calls the \command{predict}
  function.

  Note that all performance values (prediction error, variable
  importance) are standardized by the overall y-standard deviation.
  Thus, reported RMSE (root-mean-squared-error) is actually standardized
  RMSE.  Values are reported at the optimal stopping time.
}
\value{
     An object of class \code{(boostmtree, predict)}, which is a list with the
     following components:
     \item{obj}{The original boosting object.}
     \item{x}{The test x-values, but with only one row per individual
       (i.e. duplicated rows are removed).}
     \item{time}{List with each component containing the time
       points for a given test individual.}
     \item{time.unq}{Vector of the unique time points.}
     \item{y}{List containing the test y-values.}
     \item{beta}{Estimated beta coefficients at the optimized \code{M}.}
     \item{mu}{Estimated mean profile at the optimized \code{M}.}
     \item{muhat}{Extrapolated mean profile to all unique time points
       evaluated at the the optimized \code{M}.}
     \item{phi}{Variance parameter at the optimized \code{M}.}
     \item{rho}{Correlation parameter at the optimized \code{M}.}
     \item{err.rate}{Test set standardized l1-error and RMSE.}
     \item{rmse}{Test set standardized RMSE at the optimized \code{M}.}
     \item{vimp}{RMSE standardized variable importance (VIMP) values.}
     \item{Mopt}{The optimized \code{M}.}
}
\author{
  Hemant Ishwaran, Pande A. and Udaya B. Kogalur
}
\references{
  Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
  Blackstone E.H., Ishwaran H. (2017).  Boosted multivariate trees for
  longitudinal data, \emph{Machine Learning}, 106(2): 277--305. 
}
\seealso{
  \command{\link{plot.boostmtree}},
  \command{\link{print.boostmtree}}
}
\examples{
\dontrun{
##------------------------------------------------------------
## Synthetic example
##
##  high correlation, quadratic time with quadratic interaction
##  largish number of noisy variables
##
##  illustrates how modified gradient improves performance
##  also compares performance to ideal and well specified linear models 
##----------------------------------------------------------------------------

## simulate the data
## simulation 2: main effects (x1, x3, x4), quad-time-interaction (x2)
dtaO <- simLong(n = 100, ntest = 100, model = 2, q = 25)

## save the data as both a list and data frame
dtaL <- dtaO$dtaL
dta <- dtaO$dta

## get the training data
trn <- dtaO$trn

## save formulas for linear model comparisons
f.true <- dtaO$f.true
f.linr <- "y~x1+x2+x3+x4+x1*time+x2*time+x3*time+x4*time"


## modified tree gradient (default)
o.1 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn],dtaL$y[trn],
         M = 350)
p.1 <- predict(o.1, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])

## non-modified tree gradient (nmtg)
o.2 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
             M = 350, mod.grad = FALSE)
p.2 <- predict(o.2, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])

## set rho = 0
o.3 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
             M = 350, rho = 0)
p.3 <- predict(o.3, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])


##rmse values compared to generalized least squares (GLS)
##for true model and well specified linear models (LM)
cat("true LM           :", boostmtree:::gls.rmse(f.true,dta,trn),"\n")
cat("well specified LM :", boostmtree:::gls.rmse(f.linr,dta,trn),"\n")
cat("boostmtree        :", p.1$rmse,"\n")
cat("boostmtree  (nmtg):", p.2$rmse,"\n")
cat("boostmtree (rho=0):", p.3$rmse,"\n")

##predicted value plots
plot(p.1)
plot(p.2)
plot(p.3)
}
}
\keyword{predict}
\keyword{boosting}

