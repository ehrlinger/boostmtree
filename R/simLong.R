####**********************************************************************
####**********************************************************************
####
####  BOOSTED MULTIVARIATE TREES FOR LONGITUDINAL DATA (BOOSTMTREE)
####  Version 1.5.1 (_PROJECT_BUILD_ID_)
####
####  Copyright 2016, University of Miami
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 3
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Project Partially Funded By:
####  ----------------------------------------------------------------
####  Dr. Ishwaran's work was funded in part by grant R01 CA163739 from
####  the National Cancer Institute.
####
####  Dr. Kogalur's work was funded in part by grant R01 CA163739 from
####  the National Cancer Institute.
####  ----------------------------------------------------------------
####  Written by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Amol Pande, Ph.D.
####    Assistant Staff,
####    Thoracic and Cardiovascular Surgery
####    Heart and Vascular Institute
####    JJ4, Room 508B,
####    9500 Euclid Ave,
####    Cleveland Clinic, Cleveland, Ohio, 44195
####
####    email:  amoljpande@gmail.com
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************










#' Simulate longitudinal data
#' 
#' Simulates longitudinal data with continuous or binary response from models
#' with increasing complexity of covariate-time interactions.
#' 
#' Simulates longitudinal data with 3 main effects and (possibly) a
#' covariate-time interaction.  Complexity of the model is specified using the
#' option \code{model}: \enumerate{ \item \emph{\code{model=0}:} Linear with no
#' covariate-time interactions.  \item \emph{\code{model=1}:} Linear
#' covariate-time interaction.  \item \emph{\code{model=2}:} Quadratic
#' time-quadratic covariate interaction.  \item \emph{\code{model=3}:}
#' Quadratic time-quadratic two-way covariate interaction.  } For details see
#' Pande et al. (2017).
#' 
#' @param n Requested training sample size.
#' @param ntest Requested test sample size.
#' @param N Parameter controlling number of time points per subject.
#' @param rho Correlation parameter.
#' @param type Type of correlation matrix.
#' @param model Requested simulation model.
#' @param family Family of response \code{y}. Use any one from {"Continuous",
#' "Binary"} based on the scale of \code{y}.
#' @param phi Variance of measurement error.
#' @param q Number of zero-signal variables (i.e., variables unrelated to y).
#' @param ... Further arguments passed to or from other methods.
#' @return An invisible list with the following components: \item{dtaL}{List
#' containing the simulated data in the following order: \code{features},
#' \code{time}, \code{id} and \code{y}.} \item{dta}{Simulated data given as a
#' data frame.} \item{trn}{Index of \code{id} values identifying the training
#' data.} \item{f.true}{Formula of the simulation model.}
#' @author Hemant Ishwaran, Amol Pande and Udaya B. Kogalur
#' @references Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
#' Blackstone E.H., Ishwaran H. (2017).  Boosted multivariate trees for
#' longitudinal data, \emph{Machine Learning}, 106(2): 277--305.
#' @keywords simulation variable selection
#' @examples
#' 
#' \dontrun{
#' ##------------------------------------------------------------
#' ##  Response is continuous
#' ##----------------------------------------------------------------------------
#' 
#' ## set the number of boosting iterations
#' M <- 500
#' 
#' ## simulation 0: only main effects (x1, x3, x4)
#' dta <- simLong(n = 100, ntest = 100, model = 0, family = "Continuous", q = 5)
#' trn <- dta$trn
#' dtaL <- dta$dtaL
#' dta <- dta$dta
#' obj.0 <-  boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn], 
#'           family = "Continuous", M = M)
#' pred.0 <- predict(obj.0, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])
#' 
#' 
#' 
#' ##------------------------------------------------------------
#' ##  Response is binary
#' ##----------------------------------------------------------------------------
#' 
#' ## set the number of boosting iterations
#' M <- 500
#' 
#' ## simulation 0: only main effects (x1, x3, x4)
#' dta <- simLong(n = 100, ntest = 100, model = 0, family = "Binary", q = 5)
#' trn <- dta$trn
#' dtaL <- dta$dtaL
#' dta <- dta$dta
#' obj.0 <-  boostmtree(dtaL$features[trn, ], dtaL$time[trn], dtaL$id[trn], dtaL$y[trn], 
#'           family = "Binary", M = M)
#' pred.0 <- predict(obj.0, dtaL$features[-trn, ], dtaL$time[-trn], dtaL$id[-trn], dtaL$y[-trn])
#' }
#' 
simLong <-  function(n = 100,
                     ntest = 0,
                     N = 5,
                     rho = 0.8,
                     type = c("corCompSym", "corAR1", "corSymm", "iid"),
                     model = c(0, 1, 2, 3),
                     family = c("Continuous", "Binary"),
                     phi = 1,
                     q = 0,
                     ...) {
  if (length(family) != 1) {
    stop("Specify any one of the family")
  }
  if (any(is.na(match(family, c(
    "Continuous", "Binary"
  ))))) {
    stop("family must be Continuous or Binary")
  }
  type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
  model <- as.numeric(match.arg(as.character(model), as.character(0:3)))
  dta <- data.frame(do.call("rbind", lapply(1:(n + ntest), function(i) {
    Ni <- round(runif(1, 1, 3 * N))
    type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
    if (type == "corCompSym") {
      corr.mat <- matrix(rho, nrow = Ni, ncol = Ni)
      diag(corr.mat) <- 1
    }
    if (type == "corAR1") {
      corr.mat <- diag(rep(1, Ni))
      if (Ni > 1) {
        for (ii in 1:(Ni - 1)) {
          corr.mat[ii, (ii + 1):Ni] <- rho^(1:(Ni - ii))
        }
        ind <- lower.tri(corr.mat)
        corr.mat[ind] <- t(corr.mat)[ind]
      }
    }
    if (type == "iid") {
      corr.mat <- diag(rep(1, Ni))
    }
    eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
    x1 <- rnorm(1)
    x2 <- runif(1, 1, 2)
    x3 <- runif(1, 1, 3)
    x4 <- rnorm(1)
    x <- c(x1, x2, x3, x4)
    p <- length(x)
    if (q > 0) {
      xnoise <- rnorm(q)
      x <- c(x, xnoise)
    }
    tm <- sample((1:(3 * N)) / N, size = Ni, replace = TRUE)
    if (model == 0) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .6 * x4 + eps
    }
    if (model == 1) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * tm  * x2   + eps
    }
    if (model == 2) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * (tm^2) * (x2^2) +
        eps
    }
    if (model == 3) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * exp(x4) -
        .65 * (tm^2) * (x2^2) * x3  + eps
    }
    mu <- GetMu(Linear_Predictor = l_pred, Family = family)
    if (family == "Continuous") {
      y <- mu
    } else {
      y <- rbinom(n = Ni,
                  size = 1,
                  prob = mu)
    }
    cbind(matrix(
      x,
      nrow = Ni,
      ncol = length(x),
      byrow = TRUE
    ),
    tm,
    rep(i, Ni),
    y)
  })))
  d <- q + 4
  colnames(dta) <- c(paste("x", 1:d, sep = ""), "time", "id", "y")
  dtaL <- list(
    features = dta[, 1:d],
    time = dta$time,
    id = dta$id,
    y = dta$y
  )
  if (model == 0) {
    f.true <- "y ~ g( x1 + x3 + x4 )"
  }
  if (model == 1) {
    f.true <- "y ~ g( x1 + x3 + x4 + I(time * x2) )"
  }
  if (model == 2) {
    f.true <- "y ~ g( x1 + x3 + x4 + I(time^2 * x2^2) )"
  }
  if (model == 3) {
    f.true <- "y ~ g( x1 + x3 + exp(x4) + I(time^2 * x2^2 * x3) )"
  }
  trn <- c(1:sum(dta$id <= n))
  return(invisible(list(
    dtaL = dtaL,
    dta = dta,
    trn = trn,
    f.true = f.true
  )))
}
