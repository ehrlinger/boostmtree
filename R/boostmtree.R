####**********************************************************************
####**********************************************************************
####
####  BOOSTED MULTIVARIATE TREES FOR LONGITUDINAL DATA (BOOSTMTREE)
####  Version 1.0.0 (_PROJECT_BUILD_ID_)
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
####    Amol Pande
####    Division of Biostatistics
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  amoljpande@gmail.com
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Consultant Staff
####    Deptartment of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####
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


boostmtree <- function(x,
                       tm,
                       id,
                       y,
                       M = 200,
                       nu = 0.05,
                       K = 5,
                       nknots = 10,
                       d = 3,
                       pen.ord = 3,
                       lambda,
                       lambda.max = 1e6,
                       lambda.iter = 2,
                       importance = FALSE,
                       svd.tol = 1e-6,
                       forest.tol = 1e-3,
                       verbose = TRUE,
                       ...)
{
  univariate <- FALSE
  id.unq <- sort(unique(id))
  n <- length(id.unq)
  if (length(id.unq) == length(id)) {
    univariate <- TRUE
    tm <- rep(0, n)
    d <- -1
  }
  user.option <- list(...)
  if (any(is.na(id)) || any(is.na(y)) || any(is.na(x)) || any(is.na(tm))) {
    stop("missing values encountered: remove observations with missing values")
  }
  x <- as.data.frame(x)
  p <- ncol(x)
  xvar.names <- colnames(x)
  X <- do.call(rbind, lapply(1:n, function(i) {
    x[id == id.unq[i],, drop = FALSE][1,, drop = FALSE]}))
  x <- do.call(rbind, lapply(1:n, function(i) {x[id == id.unq[i],, drop = FALSE]}))
  Ymean <- mean(y, na.rm = TRUE)
  Ysd <- sd(y, na.rm = TRUE)
  if (Ysd < 1e-6) {
    Ysd <- 1
  }
  Yorg <- lapply(1:n, function(i) {y[id == id.unq[i]]})
  Y <- lapply(1:n, function(i) {(y[id == id.unq[i]] - Ymean) / Ysd})
  ni <- unlist(lapply(1:n, function(i) {sum(id == id.unq[i])}))
  id <- sort(id)
  tm.unq <- sort(unique(tm))
  n.tm <- length(tm.unq)
  tm.id <- lapply(1:n, function(i) {
    tm.i <- tm[id == id.unq[i]]
    match(tm.i, tm.unq)
  })
  if (nknots < 0) {
    warning("bsplines require a positive number of knots: eliminating b-spline fitting")
    d <- 0
  }
  if (d >= 1) {
    if (n.tm > 1) {
      bs.tm <- bs(tm.unq, df = nknots + d, degree = d)
      X.tm <- cbind(1, bs.tm)
      attr(X.tm, "knots") <- attr(bs.tm, "knots")
      attr(X.tm, "Boundary.knots") <- attr(bs.tm, "Boundary.knots")
    }
    else {
      X.tm <- cbind(1, cbind(tm.unq))
    }
  }
  else {    
    X.tm <- cbind(rep(1, n.tm))
    lambda <- 0
  }
  df.D <- ncol(X.tm)
  D <- lapply(1:n, function(i) {
    cbind(X.tm[tm.id[[i]],, drop = FALSE])
  })
  nu <- {if (length(nu) > 1) nu else rep(nu, 2)}
  if (sum(!(0 < nu & nu <= 1)) > 0) {
    stop("regularization parameter (nu) must be in (0,1]")
  }
  nu.vec <- c(nu[1], rep(nu[2], df.D - 1))
  ntree <- is.hidden.ntree(user.option)
  bootstrap <- is.hidden.bootstrap(user.option)
  if (ntree == 1) {
    nodesize <- max(1, round(n/(2 * K)))
    mtry <- df.D + p
  }
  else {
    nodedepth <- max(0, log(max(0, K - 1), base = 2))
    nodesize <- 1
    mtry <- NULL
    if (missing(lambda) || lambda < 0) {
      lambda <- 0
    }
  }
  importance <- FALSE
  vimpFlag <- bootstrap == "by.root" && importance && ntree == 1
  vimp <- NULL
  if (ntree > 1) {
    learnerUsed <- "mforest"
  }
  else {
    if (df.D == 1) {
      learnerUsed <- "mtree"
    }
    else {
      learnerUsed <- "mtree-Pspline"
    }
  }
  lambda.est.flag <- FALSE
  pen.lsq.matx <- penBSderiv(df.D - 1, pen.ord)
  if (!univariate && ntree == 1 && (missing(lambda) || lambda < 0)) {
    if (df.D >= (pen.ord + 2)) {
      lambda.est.flag <- TRUE
      pen.mix.matx <- penBS(df.D - 1, pen.ord)
      svd.pen <- svd(pen.mix.matx)
      d.zap <- svd.pen$d < svd.tol
      d.sqrt <- sqrt(svd.pen$d)
      d.sqrt[d.zap] <- 0
      d.inv.sqrt <- 1 / sqrt(svd.pen$d)
      d.inv.sqrt[d.zap] <- 0
      pen.inv.sqrt.matx <- svd.pen$v %*% (t(svd.pen$v) * d.inv.sqrt)
    }
    else {
      warning("not enough degrees of freedom to estimate lambda: setting lambda to zero\n")
      lambda <- 0
    }
  }
  mu <- lapply(1:n, function(i) {rep(0, ni[i])})
  beta <- matrix(0, n, df.D)
  if (vimpFlag) {
    vimp <- matrix(0, M, p)
    colnames(vimp) <- xvar.names
  }
  if (ntree == 1) {
    baselearner <- membership.list <- gamma.list <- vector("list", length = M)
  }
  else {
   membership.list <- gamma.list <- NULL
   baselearner <- vector("list", length = M)
  }
  lambda.vec <- phi.vec <- rho.vec <- rep(0, M)
  lambda.initial <- var(unlist(Y), na.rm = TRUE)
  rho.fit.flag <- is.hidden.rho(user.option)
  if (rho.fit.flag == TRUE) {
    rho <- 0
  }
  else {
    rho <- rho.fit.flag
    rho.fit.flag <- FALSE
    if (rho < 0 || rho > 1) {
      stop("user specified rho is not valid:", rho)
    }
  }
  sigma <- phi <- 1
  if (!lambda.est.flag) {
    sigma <- lambda * (1 - rho)
  }
  Y.names <- paste("Y", 1:df.D, sep = "")
  rfsrc.f <- as.formula(paste("Multivar(", paste(Y.names, collapse = ","), paste(") ~ ."), sep = ""))
  if (verbose) cat("  implementing multivariate boosting...\n")
  for (m in 1:M) {
    if (verbose) {
      cat("\t-- iteration:", m, "\n")
    }
    gm <- t(matrix(unlist(lapply(1:n, function(i) {
      rmi <- rho.inv(ni[i], rho)
      cmi <- rmi * sum(Y[[i]] - mu[[i]], na.rm = TRUE)
      t(D[[i]]) %*% (Y[[i]] - mu[[i]] - cbind(rep(cmi, ni[i])))
    })), nrow = df.D))
    incoming.data <- cbind(gm, X)
    names(incoming.data) = c(Y.names, names(X))
    if (ntree > 1) {
      rfsrc.obj <- rfsrc(rfsrc.f,
                         data = incoming.data,
                         mtry = mtry,
                         nodedepth = nodedepth,
                         nodesize = nodesize,
                         importance = "none",
                         bootstrap = bootstrap,
                         ntree = ntree,
                         forest.wt = TRUE)
      Kmax <- max(rfsrc.obj$leaf.count, na.rm = TRUE)
      baselearner[[m]] <- list(forest = rfsrc.obj)
    }
    else {
      rfsrc.obj <- rfsrc(rfsrc.f,
                         data = incoming.data,
                         ntree = 1,
                         mtry = mtry,
                         nodesize = nodesize,
                         importance = "none",
                         bootstrap = bootstrap)
      baselearner[[m]] <- rfsrc.obj
      result.pred <- predict.rfsrc(rfsrc.obj,
                                   membership = TRUE,
                                   ptn.count = K,
                                   importance = "none")
      membership <- membership.org <- c(result.pred$ptn.membership)
      membership.list[[m]] <- membership.org
      membership <- as.numeric(factor(membership))
      ptn.id <- unique(membership)
      Kmax <-  length(ptn.id)
      if (vimpFlag) {
        oob <- which(rfsrc.obj$inbag == 0)
        n.oob <- length(oob)
        Xnoise <- do.call(rbind, lapply(1:p, function(k) {
          X.k <- X[oob,, drop = FALSE]
          X.k[, k] <- sample(X.k[, k])
          X.k
        }))
        membershipNoise <- c(predict.rfsrc(rfsrc.obj,
                        newdata = Xnoise,
                        membership = TRUE,
                        ptn.count = K,
                        importance = "none")$ptn.membership)
        membershipNoise <- as.numeric(factor(membershipNoise, levels = levels(factor(membership.org))))
      }
    }
    if (ntree == 1) {
      if (lambda.est.flag) {
        transf.data <- mclapply(1:n, function(i) {
          if (ni[i] > 1) {
            ci <- rho.inv.sqrt(ni[i], rho)
            R.inv.sqrt <- (diag(1, ni[i]) - matrix(ci, ni[i], ni[i])) / sqrt(1 - rho)
          }
          else {
            R.inv.sqrt <- cbind(1)
          }
          Ynew <- R.inv.sqrt %*% (Y[[i]] - mu[[i]])
          Xnew <- R.inv.sqrt %*% D[[i]][, 1, drop = FALSE]
          Znew <- R.inv.sqrt %*% D[[i]][, -1, drop = FALSE] %*% pen.inv.sqrt.matx
          list(Ynew = Ynew, Xnew = Xnew, Znew = Znew)
        })
        lambda.hat <- lambda.initial
        for (k in 1:lambda.iter) {
          blup.obj <-  blup.solve(transf.data, membership, lambda.hat, Kmax)
          lambda.obj <- lapply(1:Kmax, function(k) {
            pt.k <- (membership == k)
            Z <- do.call(rbind, lapply(which(pt.k), function(j) {transf.data[[j]]$Znew}))
            X <- do.call(rbind, lapply(which(pt.k), function(j) {transf.data[[j]]$Xnew}))
            Y <- unlist(lapply(which(pt.k), function(j) {transf.data[[j]]$Ynew}))
            ZZ <- t(Z) %*% Z
            rss <- (Y - X %*% c(blup.obj[[k]]$fix.eff))^2
            robust.pt <- (rss <= quantile(rss, .99, na.rm = TRUE))
            rss <- sum(rss[robust.pt], na.rm = TRUE)
            resid <- (Y - X %*% c(blup.obj[[k]]$fix.eff) - Z %*% c(blup.obj[[k]]$rnd.eff))^2
            resid <- resid[robust.pt]
            return(list(trace.Z = sum(diag(ZZ)), rss = rss, resid = resid))
          })
          num <- sum(unlist(lapply(1:Kmax, function(k) {lambda.obj[[k]]$trace.Z})), na.rm = TRUE)
          den <- sum(unlist(lapply(1:Kmax, function(k) {lambda.obj[[k]]$rss})), na.rm = TRUE)
          N <- sum(unlist(lapply(1:Kmax, function(k) {lambda.obj[[k]]$resid})), na.rm = TRUE)
          if (!is.na(den) && den > (.99 * N)) {
            lambda.hat <- num / (den - .99 * N)
          }
          else {
            lambda.hat <- min(lambda.hat, lambda.max)
          }
          lambda.hat <- min(lambda.hat, lambda.max)
        }
        lambda <- lambda.hat 
        sigma <- lambda * (1 - rho) 
      }
      Xnew <- mclapply(1:n, function(i) {
        rmi <- rho.inv(ni[i], rho)
        Wi <- diag(1, ni[i]) - matrix(rmi, ni[i], ni[i])
        t(D[[i]]) %*% Wi %*% D[[i]]
      })
      gamma <- lapply(1:Kmax, function(k) {
        pt.k <- (membership == k)
        if (sum(pt.k) > 0) {
          YnewSum <- colSums(gm[pt.k,, drop = FALSE])
          XnewSum <- Reduce("+", lapply(which(pt.k), function(j) {Xnew[[j]]}))
          XnewSum <- XnewSum + sigma * pen.lsq.matx
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
      })
      gamma.matx <- matrix(0, Kmax, df.D+1)
      gamma.matx[, 1] <- sort(unique(membership.org))
      gamma.matx[, 2:(df.D+1)] <- matrix(unlist(gamma), ncol = df.D, byrow = TRUE)
      gamma.list[[m]] <- gamma.matx
      beta.update <- matrix(unlist(lapply(1:n, function(i) {
        gamma[[membership[i]]]})), nrow = df.D)
      beta.old <- beta
      beta <- beta.old + t(beta.update * nu.vec)
      mu.old <- mu
      mu <- lapply(1:n, function(i) {D[[i]] %*% beta[i, ]})
      if (vimpFlag) {
        beta.std <- beta[oob,, drop = FALSE] * Ysd
        beta.std[, 1] <- beta.std[, 1] + Ymean
        err.oob <- l2Dist(Yorg[oob], lapply(1:n.oob, function(i) {D[[oob[i]]] %*% beta.std[i, ]}))
        beta.update.vimp <- lapply(1:p, function(k) {
          membership.k <- membershipNoise[((k-1) * n.oob + 1):(k * n.oob)]
          matrix(unlist(lapply(1:n.oob, function(i) {
            gamma[[membership.k[i]]]})), nrow = df.D)
        })
        err.vimp <- sapply(1:p, function(k) {
          beta.vimp.k <- beta.old[oob,, drop = FALSE] + t(beta.update.vimp[[k]] * nu.vec)
          beta.vimp.k <- beta.vimp.k * Ysd
          beta.vimp.k[, 1] <- beta.vimp.k[, 1] + Ymean
          l2Dist(Yorg[oob], lapply(1:n.oob, function(i) {D[[oob[i]]] %*% beta.vimp.k[i, ]}))
        })
        vimp[m, ] <- 100 * (err.vimp - err.oob) / Ysd
      }
    }
    else{
      forest.wt <- rfsrc.obj$forest.wt
      Xnew <- mclapply(1:n, function(i) {
        rmi <- rho.inv(ni[i], rho)
        Wi <- diag(1, ni[i]) - matrix(rmi, ni[i], ni[i])
        t(D[[i]]) %*% Wi %*% D[[i]]
      })
      bhat <- do.call("cbind", mclapply(1:n, function(i) {
        fwt.i <- forest.wt[i, ]
        fwt.i[fwt.i <= forest.tol] <- 0
        pt.i <- (fwt.i != 0)
        if (sum(pt.i) > 0) {
          fwt.i <- fwt.i / sum(fwt.i)
          YnewSum <- colSums(fwt.i[pt.i] * gm[pt.i,, drop = FALSE])
          XnewSum <- Reduce("+", lapply(which(pt.i), function(j) {fwt.i[j] * Xnew[[j]]}))
          XnewSum <- XnewSum + sigma * pen.lsq.matx
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
      beta <- beta + t(bhat * nu.vec)
      mu <- lapply(1:n, function(i) {D[[i]] %*% beta[i, ]})
      baselearner[[m]] <- c(baselearner[[m]],
                            list(gm = gm),
                            list(Xnew = Xnew),
                            list(pen = sigma * pen.lsq.matx))
    }
    if (!univariate) {
      resid.data <- data.frame(y  = unlist(lapply(1:n, function(i) {Y[[i]] - mu[[i]]})),
                               x,
                               tm = unlist(lapply(1:n, function(i) {tm[id == id.unq[i]]})),
                               id = unlist(lapply(1:n, function(i) {rep(id.unq[i], ni[i])})))
      gls.obj <- tryCatch({gls(y ~ ., data = resid.data,
                               correlation = corCompSymm(form = ~ 1 | id))},
                          error = function(ex){NULL})
      if (is.null(gls.obj)) {
        gls.obj <- tryCatch({gls(y ~ 1, data = resid.data,
                                 correlation = corCompSymm(form = ~ 1 | id))},
                            error = function(ex){NULL})
      }
      if (!is.null(gls.obj)) {
        phi <- gls.obj$sigma^2
        if (rho.fit.flag) {
          rho <- as.numeric(coef(gls.obj$modelStruct$corStruc, unconstrained = FALSE))
          rho <- max(min(0.999, rho, na.rm = TRUE), -0.999)
        }
      }
      phi.vec[m] <- phi * Ysd ^ 2
      rho.vec[m] <- rho
      if (verbose) {
        cat("phi   :", phi.vec[m], "\n")
        cat("rho   :", rho.vec[m], "\n")
      }
      sigma <- lambda * (1 - rho) 
      lambda.vec[m] <- lambda
      if (verbose) {
        cat("lambda:", lambda.vec[m], "\n")
      }
    }
  }
  beta <- beta * Ysd
  beta[, 1] <- beta[, 1] + Ymean
  mu <- lapply(1:n, function(i) {c(mu[[i]] * Ysd + Ymean)})
  y <- lapply(1:n, function(i) {y[id == id.unq[i]]})
  obj <- list(x = X,
              xvar.names = xvar.names,
              time = lapply(1:n, function(i) {tm[id == id.unq[i]]}),
              id = id,
              y = Yorg,
              ymean = Ymean,
              ysd = Ysd,
              gamma = gamma.list,
              beta = beta,
              mu = mu,
              lambda = lambda.vec,
              phi = phi.vec,
              rho = rho.vec,
              baselearner = baselearner,
              membership = membership.list,
              vimp = if (!is.null(vimp)) colSums(vimp, na.rm = TRUE) else NULL,
              D = X.tm,
              d = d,
              pen.ord = pen.ord,
              K = K,
              M = M,
              nu = nu,
              ntree = ntree)
  class(obj) <- c("boostmtree", "grow", learnerUsed)
  invisible(obj)
}
