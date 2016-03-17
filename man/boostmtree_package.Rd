\name{boostmtree-package}
\alias{boostmtree-package}
\docType{package}
\title{
  Boosted multivariate trees for longitudinal data.
}
\description{
  Multivariate extension of Friedman's (2001) gradient descent boosting
  method to longitudinal data using multivariate tree base learners.
  Covariate-time interactions are modeled using penalized B-splines
  (P-splines) with estimated adaptive smoothing parameter.
}

\section{Package Overview}{
  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.

  \enumerate{
    \item \command{\link{boostmtree}}

    This is the main entry point to the package.  It grows a
    multivariate tree using user supplied training data.  Trees are
    grown using the \pkg{randomForestSRC} R-package.

    \item \command{\link{predict.boostmtree}} (\command{predict})

    Used for prediction.  Predicted values are obtained by dropping the
    user supplied test data down the grow forest.  The resulting object
    has class \code{(rfsrc, predict)}.

  }
}
\author{
  Hemant Ishwaran, Pande A. and Udaya B. Kogalur
}
\references{
    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Friedman J.H. (2001). Greedy function approximation: a gradient
    boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.
     
    Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.

    Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.

    Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2014). Random survival forests for competing risks.  To appear
    in \emph{Biostatistics}.

    Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
    Blackstone E.H., Ishwaran H. (2016).  Boosted multivariate trees for
    longitudinal data.
}
\seealso{
  \command{\link{partialPlot}},
  \command{\link{plot.boostmtree}},
  \command{\link{predict.boostmtree}},
  \command{\link{print.boostmtree}},
  \command{\link{simLong}}
}
\keyword{package}


