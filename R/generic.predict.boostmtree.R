####**********************************************************************
####**********************************************************************
####
####  BOOSTED MULTIVARIATE TREES FOR LONGITUDINAL DATA (BOOSTMTREE)
####  Version 1.4.0 (_PROJECT_BUILD_ID_)
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


generic.predict.boostmtree <- function(object,
                                       x,
                                       tm,
                                       id,
                                       y,
                                       M,
                                       eps = 1e-5,
                                       ...)
{
  if (missing(object)) {
    stop("object is missing!")
  }
  if (sum(inherits(object, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }
  user.option <- match.call(expand.dots = TRUE)
  partial <- is.hidden.partial(user.option)
  if (!partial) {
    if (missing(x)) {
      X <- object$x
      n <- nrow(X)
      X.tm <- object$X.tm
      tm <- object$time
      tm.unq <- sort(unique(unlist(object$time)))
      id <- object$id
      id.unq <- sort(unique(id))
      testFlag <- FALSE
    }
    else {
      if (!missing(x) && (missing(id) || missing(tm))) {
        X <- x
        n <- nrow(X)
        tm.unq <- sort(unique(unlist(object$time)))
        tm <- lapply(1:n, function(i){tm.unq})
        id <- id.unq <- 1:n
        if (!missing(y)) {
          Y <- lapply(1:n, function(i) {y[i]})
          testFlag <- TRUE
        }
        else {
          testFlag <- FALSE
        }
      }
      else{
        if (missing(id)) {
          stop("test set id values are missing\n")
        }
        id.unq <- sort(unique(id))
        n <- length(id.unq)
        if (missing(x)) {
          stop("test set x values are missing\n")
        }
        X <- do.call(rbind, lapply(1:n, function(i) {
          x[id == id.unq[i],, drop = FALSE][1,, drop = FALSE]}))
        if (missing(tm)) {
          stop("test set time values are missing\n")
        }
        tm.unq <- sort(unique(tm))
        if (!missing(y)) {
          tm <- lapply(1:n, function(i) {tm[id == id.unq[i]]})
          Y <- lapply(1:n, function(i) {y[id == id.unq[i]]})
          testFlag <- TRUE
        }
        else {
          testFlag <- FALSE
        }
      }
    }
  }
  else {
      X <- x
      n <- nrow(X)
      tm.unq <- tm
      tm <- lapply(1:n, function(i){tm.unq})
      id.unq <- 1:n
      id <- rep(id.unq,each = length(tm.unq))
      testFlag <- FALSE
  }
  if (missing(M)) {
    M <- object$M
    Mflag <- FALSE
  }
  else {
    M <- max(1, min(M, object$M))
    Mflag <- TRUE
  }
  K <- object$K
  nu <- object$nu
  ntree <- object$ntree
  p <- ncol(X)
  X.tm <- object$X.tm
  df.D <- ncol(X.tm)
  xvar.names <- colnames(X)
  nu.vec <- c(nu[1], rep(nu[2], df.D - 1))
  Ymean <- object$ymean
  Ysd <- object$ysd
  gamma <- object$gamma
  baselearner <- object$baselearner
  family <- object$family
  id.index <- lapply(1:n,function(i){
    which(id == id.unq[i] )
  })
  ni <- unlist(lapply(1:n,function(i){
    length(id.index[[i]])
  }))
  tm.unq_train <- object$tm.unq
  tm.index <- lapply(1:n,function(i){
    AppoxMatch(unlist(tm)[id.index[[i]]],tm.unq_train)
  })
  D <- lapply(1:n,function(i){
    X.tm[tm.index[[i]] ,,drop = FALSE]
  })
  beta <- matrix(0, n, df.D)
  if (ntree == 1) {
    beta.vimp <- beta.cov.vimp <- beta.time.vimp <- vector("list", p)
  }
  l_pred.list <- vector("list", M)
  forest.tol <- object$forest.tol
  if (ntree == 1) {
    rf.cores.old <- getOption("rf.cores")
    mc.cores.old <- getOption("mc.cores")
    membership <- papply(1:M, function(m) {
      options(rf.cores = 1, mc.cores = 1)
      c(predict.rfsrc(baselearner[[m]],
                      newdata = X,
                      membership = TRUE,
                      ptn.count = K,
                      na.action = "na.impute",
                      importance = "none")$ptn.membership)
    })
    nullObj <- lapply(1:M, function(m) {
      orgMembership <- gamma[[m]][, 1]
      if (m == 1) {
        l_pred.list[[m]] <<- lapply(1:n, function(i) {
          D[[i]]%*%t(gamma[[m]][match(membership[[m]][i], orgMembership), -1, drop = FALSE]*nu.vec) 
        })
      }
      else {
        l_pred.list[[m]] <<- lapply(1:n, function(i) {
          unlist(l_pred.list[[m-1]][i]) + D[[i]]%*%t(gamma[[m]][match(membership[[m]][i], orgMembership), -1, drop = FALSE]*nu.vec)
        })
      }
      NULL
    })
    rm(nullObj)
    l_pred.list <- lapply(l_pred.list, function(mlist){  
      lapply(1:n,function(i) {mlist[[i]] * Ysd + Ymean})
    })
    mu.list <- lapply(1:M,function(m){
      GetMu(Linear_Predictor = l_pred.list[[m]],Family = family)
    })
    if (testFlag) {
      err.rate <- matrix(unlist(lapply(mu.list, function(mlist) {
        c(l1Dist(Y, mlist), l2Dist(Y, mlist)) 
      })), ncol = 2, byrow = TRUE)
      colnames(err.rate) <- c("l1", "l2")
    }
    else {
      err.rate <- NULL
    }
    if (!Mflag && testFlag) {
      diff.err <- abs(err.rate[, "l2"] - min(err.rate[, "l2"], na.rm = TRUE))
      diff.err[is.na(diff.err)] <- 1
      if (sum(diff.err < Ysd * eps) > 0) {
        Mopt <- min(which(diff.err < eps))
      }
      else {
        Mopt <- M
      }
    }
    else {
      Mopt <- M
    }
  }
  else{
    nullObj <- lapply(1:M, function(m) {
      gm <- baselearner[[m]]$gm
      Xnew <- baselearner[[m]]$Xnew
      pen <- baselearner[[m]]$pen
      forest.wt <- predict.rfsrc(baselearner[[m]]$forest,
                                 newdata = X,
                                 importance = "none",
                                 forest.wt = TRUE)$forest.wt
      beta.m.org <- do.call("cbind", mclapply(1:n, function(i) {
        fwt.i <- forest.wt[i, ]
        fwt.i[fwt.i <= forest.tol] <- 0
        pt.i <- (fwt.i != 0)
        if (sum(pt.i) > 0) {
          fwt.i <- fwt.i / sum(fwt.i)
          YnewSum <- colSums(fwt.i[pt.i] * gm[pt.i,, drop = FALSE])
          XnewSum <- Reduce("+", lapply(which(pt.i), function(j) {fwt.i[j] * Xnew[[j]]}))
          XnewSum <- XnewSum + pen
          qr.obj <- tryCatch({qr.solve(XnewSum, YnewSum)}, error = function(ex){NULL})
          if (!is.null(qr.obj)) {
            qr.obj
          }
          else {
            rep(0, df.D)
          }
        }
        else {
          rep(0, df.D)
        }
      }))
      beta.m <- t(beta.m.org * nu.vec * Ysd)
      if (m == 1) {
        beta.m[, 1] <- beta.m[, 1] + Ymean
        beta <<- beta.m
      }
      else {
        beta <<- beta + beta.m 
      }
      Dbeta.m <- D %*% (beta.m.org * nu.vec)
      if (m == 1) {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      else {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          unlist(mu.list[[m-1]][i]) + Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      NULL
    })
    mu.list <- lapply(mu.list, function(mlist){  
      lapply(1:n,function(i) {mlist[[i]] * Ysd + Ymean})
    })
    if (testFlag) {
      err.rate <- matrix(unlist(lapply(mu.list, function(mlist) {
        c(l1Dist(Y, mlist), l2Dist(Y, mlist))
      })), ncol = 2, byrow = TRUE) 
      colnames(err.rate) <- c("l1", "l2")
    }
    else {
      err.rate <- NULL
    }
    if (!Mflag && testFlag) {
      diff.err <- abs(err.rate[, "l2"] - min(err.rate[, "l2"], na.rm = TRUE))
      diff.err[is.na(diff.err)] <- 1
      if (sum(diff.err < Ysd * eps) > 0) {
        Mopt <- min(which(diff.err < eps))
      }
      else {
        Mopt <- M
      }
    }
    else {
      Mopt <- M
    }
  }
  l_pred_hat <- lapply(1:n,function(i){
    l_pred.i <- rep(0,nrow(X.tm))
    NullObj <- lapply(1:Mopt,function(m){
        orgMembership <- gamma[[m]][, 1]
        l_pred.i <<- l_pred.i + X.tm%*%t(gamma[[m]][match(membership[[m]][i], orgMembership), -1, drop = FALSE]*nu.vec)
        NULL
    })
    l_pred.i * Ysd + Ymean
  })
  muhat <- GetMu(Linear_Predictor = l_pred_hat,Family = family)
  object$baselearner <- object$membership <- object$gamma <- NULL
  pobj <- list(boost.obj = object,
               x = X,
               time = tm,
               id = id,
               y = if (testFlag) Y else NULL,
               family = family,
               ymean = Ymean,
               ysd = Ysd,
               xvar.names = xvar.names,
               K = K,
               n = n,
               ni = ni,
               nu.vec = nu.vec,
               D = D,
               df.D = df.D,
               time.unq = tm.unq,
               baselearner = baselearner,
               gamma = gamma,
               membership = membership,
               mu = mu.list[[Mopt]],
               muhat = muhat,
               phi = object$phi[Mopt],
               rho = object$rho[Mopt],
               err.rate = if (!is.null(err.rate)) err.rate / Ysd else NULL,
               rmse = if (!is.null(err.rate)) err.rate[Mopt, "l2"] / Ysd else NULL,
               Mopt = if (testFlag) Mopt else NULL)
  class(pobj) <- c("boostmtree", "predict", class(object)[3])
  invisible(pobj)
}
