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
 proximity = FALSE,
 verbose = TRUE,
 eps = 1e-4,
 forest.tol = 1e-3,
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
  \item{proximity}{Calculate proximity of the test data to the training data?}
  \item{verbose}{Should verbose output be printed?}
  \item{eps}{Tolerance value used for determing the optimal \code{M}.
    For experts only.}
  \item{forest.tol}{Tolerance used for forest weighted least squares
      solution.  Experimental and for experts only.}
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
     \item{err.rate}{Test set l1-error and RMSE.}
     \item{mse}{Test set RMSE at the optimized \code{M}.}
     \item{vimp}{Variable importance values.}
     \item{proximity}{Proximity matrix.}
     \item{Mopt}{The optimized \code{M}.}
}
\author{
  Hemant Ishwaran, Pande A. and Udaya B. Kogalur
}
\references{
  Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
  Blackstone E.H., Ishwaran H. (2016).  Boosted multivariate trees for
  longitudinal data.
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
##  noisy variables
##
##  calculate test set error under adaptive rho and mispecified rho
##  compare to GLS (generalized least-squares)
##----------------------------------------------------------------------------

## simulate the data
dtaO <- simLong(n = 100, ntest = 100, model = 2, q = 5)

## save both the list and data frame formatted data
dtaL <- dtaO$dtaL
dta <- dtaO$dta

## identify the training data
trn <- dtaO$trn

## save formulas for GLS analysis
f.true <- dtaO$f.true
f.linr <- "y~x1+x2+x3+x4+x1*time+x2*time+x3*time+x4*time"


##adaptive rho
grow.1 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
             M = 500)
pred.1 <- predict(grow.1,
   dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])

##fixed rho = 0
grow.2 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn],
             M = 500, rho = 0)
pred.2 <- predict(grow.2,
   dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])

##mse values
cat("true model        :", boostmtree:::gls.mse(f.true,dta,trn),"\n")
cat("linear GLS        :", boostmtree:::gls.mse(f.linr,dta,trn),"\n")
cat("boostmtree        :", pred.1$mse,"\n")
cat("boostmtree (rho=0):", pred.2$mse,"\n")

##predicted value plot
plot(pred.1)
}
}
\keyword{predict}
\keyword{boosting}

