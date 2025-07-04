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











#' Print Summary Output
#' 
#' Print summary output from the boosting analysis.
#' 
#' 
#' @param x An object of class \code{(boostmtree, grow)} or \code{(boostmtree,
#' predict)}.
#' @param ... Further arguments passed to or from other methods.
#' @author Hemant Ishwaran, Amol Pande and Udaya B. Kogalur
#' @references Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
#' Blackstone E.H., Ishwaran H. (2017).  Boosted multivariate trees for
#' longitudinal data, \emph{Machine Learning}, 106(2): 277--305.
#' @keywords print
print.boostmtree <- function(x, ...) {
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &&
      sum(inherits(x, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop(
      "this function only works for objects of class `(boostmtree, grow)' or "+"'(boostmtree, predict)'"
    )
  }
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {
    univariate <- length(x$id) == length(unique(x$id))
    cat("model                       :", class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    cat("Family                      :", x$family, "\n")
    n_levels <- (x$n.Q + 1)
    if (x$family == "Nominal" || x$family == "Ordinal") {
      cat("No of levels                :", n_levels, "\n")
    }
    if (x$ntree > 1) {
      cat("ntree                     :", x$ntree, "\n")
    }
    cat("number of K-terminal nodes  :", x$K, "\n")
    cat("regularization parameter    :", x$nu[1], "\n")
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    if (!univariate) {
      cat("number of unique time points:", length(sort(unique(
        unlist(x$time)
      ))), "\n")
      cat("avg. number of time points  :", round(mean(sapply(
        x$time, length
      ), na.rm = TRUE), 2), "\n")
      cat("B-spline dimension          :", ncol(x$X.tm), "\n")
      cat("penalization order          :", x$pen.ord, "\n")
    } else {
      cat("univariate family           :", TRUE, "\n")
    }
    cat("boosting iterations         :", x$M, "\n")
    if (!is.null(x$err.rate)) {
      if (x$family == "Nominal" || x$family == "Ordinal") {
        n.Q <- x$n.Q
      } else {
        n.Q <- 1
      }
      optimized_rho <- unlist(lapply(1:n.Q, function(q) {
        if (x$family == "Nominal" || x$family == "Ordinal") {
          x$rho[x$Mopt[q], q]
        } else {
          x$rho[x$Mopt[q]]
        }
      }))
      optimized_phi <- unlist(lapply(1:n.Q, function(q) {
        if (x$family == "Nominal" || x$family == "Ordinal") {
          x$phi[x$Mopt[q], q]
        } else {
          x$phi[x$Mopt[q]]
        }
      }))
      cat("optimized number iterations :", x$Mopt, "\n")
      if (!univariate) {
        cat("optimized rho               :", round(optimized_rho, 4), "\n")
        cat("optimized phi               :", round(optimized_phi, 4), "\n")
      }
      cat("OOB cv RMSE                 :", round(x$rmse, 4), "\n")
    }
  } else {
    univariate <- length(x$boost.obj$id) == length(unique(x$boost.obj$id))
    cat("model                       :", class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    cat("Family                      :", x$family, "\n")
    n_levels <- (x$n.Q + 1)
    if (x$family == "Nominal" || x$family == "Ordinal") {
      cat("No of levels                :", n_levels, "\n")
    }
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    if (!univariate) {
      cat("number of unique time points:", length(sort(unique(
        unlist(x$time)
      ))), "\n")
      cat("avg. number of time points  :", round(mean(sapply(
        x$time, length
      ), na.rm = TRUE), 2), "\n")
      if (!is.null(x$err.rate)) {
        if (x$family == "Nominal" || x$family == "Ordinal") {
          n.Q <- x$n.Q
        } else {
          n.Q <- 1
        }
        optimized_rho <- unlist(lapply(1:n.Q, function(q) {
          if (x$family == "Nominal" || x$family == "Ordinal") {
            x$boost.obj$rho[x$Mopt[q], q]
          } else {
            x$boost.obj$rho[x$Mopt[q]]
          }
        }))
        optimized_phi <- unlist(lapply(1:n.Q, function(q) {
          if (x$family == "Nominal" || x$family == "Ordinal") {
            x$boost.obj$phi[x$Mopt[q], q]
          } else {
            x$boost.obj$phi[x$Mopt[q]]
          }
        }))
        cat("optimized number iterations :", x$Mopt, "\n")
        cat("optimized rho               :", round(optimized_rho, 4), "\n")
        cat("optimized phi               :", round(optimized_phi, 4), "\n")
        cat("test set RMSE               :", round(x$rmse, 4), "\n")
      }
    } else {
      if (!is.null(x$err.rate)) {
        cat("optimized number iterations :", x$Mopt, "\n")
        cat("test set RMSE               :", round(x$rmse, 4), "\n")
      }
      cat("univariate family           :", TRUE, "\n")
    }
  }
}
