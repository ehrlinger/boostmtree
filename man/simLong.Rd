\name{simLong}
\alias{simLong}
\title{Simulate longitudinal data}
\description{
  Simulates longitudinal data from models with increasing complexity of
  covariate-time interactions.
}
\usage{simLong(n,
 ntest = 0,
 N = 5,
 rho = 0.8,
 type = c("corCompSym", "corAR1", "corSymm", "iid"),
 model = c(0, 1, 2, 3),
 phi = 1,
 q = 0,
 ...)}
\arguments{
  \item{n}{Requested training sample size.}
  \item{ntest}{Requested test sample size.}
  \item{N}{Parameter controlling number of time points per subject.}
  \item{rho}{Correlation parameter.}
  \item{type}{Type of correlation matrix.}
  \item{model}{Requested simulation model.}
  \item{phi}{Variance of measurement error.}
  \item{q}{Number of zero-signal variables (i.e., variables unrelated to y).}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Simulates longitudinal data with 3 main effects and (possibly) a covariate-time
  interaction.  Complexity of the model is specified using the option \code{model}:
  \enumerate{
    \item \emph{\code{model=0}:} Linear with no covariate-time interactions.
    \item \emph{\code{model=1}:} Linear covariate-time interaction.
    \item \emph{\code{model=2}:} Quadratic time-quadratic covariate
    interaction.
    \item \emph{\code{model=3}:} Quadratic time-quadratic two-way
    covariate interaction.
  }
  For details see Pande et al. (2016).
}
\value{
  An invisible list with the following components:
  \item{dtaL}{List containing the simulated data in the following order:
    \code{features}, \code{time}, \code{id} and \code{y}.}
  \item{dta}{Simulated data given as a data frame.}
  \item{trn}{Index of \code{id} values identifying the training data.}
  \item{f.true}{Formula of the simulation model.}
}
\author{
  Hemant Ishwaran, Pande A. and Udaya B. Kogalur
}
\references{
  Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
  Blackstone E.H., Ishwaran H. (2016).  Boosted multivariate trees for
  longitudinal data.
}
\examples{
\dontrun{
##------------------------------------------------------------
##
##  test how well we can select variables
##----------------------------------------------------------------------------


## set the number of boosting iterations
M <- 500

## simulation 0: only main effects (x1, x3, x4)
dta <- simLong(n = 100, ntest = 100, model = 0, q = 5)
trn <- dta$trn
dtaL <- dta$dtaL
dta <- dta$dta
obj.0 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn], M = M)
pred.0 <- predict(obj.0, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])
vimp.0 <- pred.0$vimp

## simulation 3: main effects (x1, x3, x4), time-interactions (x2, x3)
dta <- simLong(n = 100, ntest = 100, model = 3, q = 5)
trn <- dta$trn
dtaL <- dta$dtaL
dta <- dta$dta
obj.3 <- boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn], M = M)
pred.3 <- predict(obj.3, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])
vimp.3 <- pred.3$vimp

## plot vimp side by side
barplot(rbind(vimp.0,vimp.3),beside=T,las=2,
    ylab="standardized vimp", args.legend = list(x = "topleft"),
    col = c(2, 4), legend.text = c("model 0", "model 3"))

}
}
\keyword{simulation}
\keyword{variable selection}
