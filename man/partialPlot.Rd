\name{partialPlot}
\alias{partialPlot}
\title{Partial plot analysis}
\description{
  Partial dependence plot of x against adjusted predicted y.
}
\usage{partialPlot(obj,
 xvar.names,
 tm,
 npts = 25,
 subset,
 conditional.xvars = NULL,
 conditional.values = NULL,
 plot.it = TRUE,
 ...)}
\arguments{
  \item{obj}{A boosting object of class \code{(boostmtree, grow)}.}
  \item{xvar.names}{Names of the x-variables to be used.  By default,
    all variables are plotted.}
  \item{tm}{Time points used for the plots of x against y.  By default,
    the deciles of the observed time values are used.}
  \item{npts}{Maximum number of points used for x.  Reduce this value if
    plots are slow.}
  \item{subset}{Vector indicating which rows of the x-data to be used for the
    analysis.  The default is to use the entire data.}
  \item{conditional.xvars}{Vector of character values indicating names of the x-variables
    to be used for further conditioning (adjusting) the predicted y values. Variable names
    should be different from \code{xvar.names}.}
  \item{conditional.values}{Vector of values taken by the variables from \code{conditional.xvars}.
    The length of the vector should be same as the length of the vector for \code{conditional.xvars},     which means only one value per conditional variable.}
  \item{plot.it}{Should plots be displayed?}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Partial dependence plot (Friedman, 2001) of x values specified by
  \code{xvar.names} against the adjusted predicted y-values over a set
  of time points specified by \code{tm}.  Analysis can be restricted to
  a subset of the data using \code{subset}. Further conditioning can be
  imposed using \code{conditional.xvars}.
}
\author{
  Hemant Ishwaran, Pande A. and Udaya B. Kogalur
}
\references{
  Friedman J.H. Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232, 2001.
}
\examples{
\dontrun{
##------------------------------------------------------------
## Synthetic example
## high correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
#simulate the data
dta <- simLong(n = 50, N = 5, rho =.80, model = 2)$dtaL

#basic boosting call
boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y, M = 300)

#plot results
#x1 has a linear main effect
#x2 is quadratic with quadratic time trend
partialPlot(boost.grow, "x1")
partialPlot(boost.grow, "x2")

#partial plot using "x2" as the conditional variable
partialPlot(boost.grow, "x1", conditional.xvar = "x2", conditional.values = 1)
partialPlot(boost.grow, "x1", conditional.xvar = "x2", conditional.values = 2)

##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------
data(spirometry, package = "boostmtree")

#boosting call: cubic B-splines with 15 knots
spr.obj <- boostmtree(spirometry$features, spirometry$time, spirometry$id, spirometry$y,
            M = 300, nu = .025, nknots = 15)

#partial plot of double-lung group at 5 years
dltx <- partialPlot(spr.obj, "AGE", tm = 5, subset=spr.obj$x$DOUBLE==1)

#partial plot of single-lung group at 5 years
sltx <- partialPlot(spr.obj, "AGE", tm = 5, subset=spr.obj$x$DOUBLE==0)

#combine the two plots: we use lowess smoothed values
dltx <- dltx$l.obj[[1]]
sltx <- sltx$l.obj[[1]]
plot(range(c(dltx[, 1], sltx[, 1])), range(c(dltx[, -1], sltx[, -1])),
     xlab = "age", ylab = "predicted y (adjusted)", type = "n")
lines(dltx[, 1], dltx[, -1], lty = 1, lwd = 2, col = "red")
lines(sltx[, 1], sltx[, -1], lty = 1, lwd = 2, col = "blue")
legend("topright", legend = c("DLTx", "SLTx"), lty = 1, fill = c(2,4))
}
}
\keyword{plot}
