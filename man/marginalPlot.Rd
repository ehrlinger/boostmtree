\name{marginalPlot}
\alias{marginalPlot}
\title{Marginal plot analysis}
\description{
  Marginal plot of x against the unadjusted predicted y. This is mainly used to obtain
  marginal relationships between x and the unadjusted predicted
  y. Marginal plots have a faster
  execution compared to partial plots (Friedman, 2001).
}
\usage{marginalPlot(obj,
 xvar.names,
 tm,
 subset,
 plot.it = TRUE,
 ...)}
\arguments{
  \item{obj}{A boosting object of class \code{(boostmtree, grow)}.}
  \item{xvar.names}{Names of the x-variables to be used.  By default,
    all variables are plotted.}
  \item{tm}{Time points used for the plots of x against y.  By default,
    the deciles of the observed time values are used.}
  \item{subset}{Vector indicating which rows of the x-data to be used for the
    analysis.  The default is to use the entire data.}
  \item{plot.it}{Should plots be displayed? If \code{xvar.names} is a vector with more than
  one variable name, then instead of displaying, plot is stored as "MarginalPlot.pdf" in the current working directory.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Marginal plot of x values specified by
  \code{xvar.names} against the unadjusted predicted y-values over a set
  of time points specified by \code{tm}.  Analysis can be restricted to
  a subset of the data using \code{subset}.
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
marginalPlot(boost.grow, "x1")
marginalPlot(boost.grow, "x2")

#Plot of all covariates. The plot will be stored as the "MarginalPlot.pdf"
# in the current working directory.
marginalPlot(boost.grow)

##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------
data(spirometry, package = "boostmtree")

#boosting call: cubic B-splines with 15 knots
spr.obj <- boostmtree(spirometry$features, spirometry$time, spirometry$id, spirometry$y,
            M = 300, nu = .025, nknots = 15)

#marginal plot of double-lung group at 5 years
dltx <- marginalPlot(spr.obj, "AGE", tm = 5, subset = spr.obj$x$DOUBLE==1)

#marginal plot of single-lung group at 5 years
sltx <- marginalPlot(spr.obj, "AGE", tm = 5, subset = spr.obj$x$DOUBLE==0)

#combine the two plots
dltx <- dltx[[1]]
sltx <- sltx[[1]]
plot(range(c(dltx[[1]][, 1], sltx[[1]][, 1])), range(c(dltx[[1]][, -1], sltx[[1]][, -1])),
     xlab = "age", ylab = "predicted y", type = "n")
lines(dltx[[1]][, 1], dltx[[1]][, -1], lty = 1, lwd = 2, col = "red")
lines(sltx[[1]][, 1], sltx[[1]][, -1], lty = 1, lwd = 2, col = "blue")
legend("topright", legend = c("DLTx", "SLTx"), lty = 1, fill = c(2,4))
}
}
\keyword{plot}
