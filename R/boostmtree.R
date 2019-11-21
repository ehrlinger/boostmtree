####**********************************************************************
####**********************************************************************
####
####  BOOSTED MULTIVARIATE TREES FOR LONGITUDINAL DATA (BOOSTMTREE)
####  Version 1.4.1 (_PROJECT_BUILD_ID_)
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


boostmtree <- function(x,
                       tm,
                       id,
                       y,
                       family = c("Continuous","Binary"),
                       M = 200,
                       nu = 0.05,
                       K = 5,
                       nknots = 10,
                       d = 3,
                       pen.ord = 3,
                       lambda,
                       lambda.max = 1e6,
                       lambda.iter = 2,
                       svd.tol = 1e-6,
                       forest.tol = 1e-3,
                       verbose = TRUE,
                       cv.flag = FALSE,
                       eps = 1e-5,
                       mod.grad = TRUE,
                       NR.iter = 3,
                       ...)
{
  if(Sys.info()["sysname"] == "Windows")
  {
    options(rf.cores=1, mc.cores= 1)
  } 
  if(length(family) != 1){
    stop("Specify any one of the family")
  }
  if(any(is.na( match(family,c("Continuous","Binary"))))){
    stop("family must be Continuous or Binary")
  }
  univariate <- FALSE
  if (missing(tm)) {
    id <- 1:nrow(x)
  }
  id.unq <- sort(unique(id))
  n <- length(id.unq)
  if (length(id.unq) == length(id)) {
    univariate <- TRUE
    tm <- rep(0, n)
    d <- -1
  }
  if (univariate) {
    mod.grad <- FALSE
    rho <- 0
    lambda.vec <- phi.vec <- rho.vec <- NULL
  }
  user.option <- list(...)
  if (any(is.na(id)) || any(is.na(y)) || any(is.na(tm))) {
    stop("missing values encountered y or id or tm: remove observations with missing values")
  }
  x <- as.data.frame(x)
  X <- do.call(rbind, lapply(1:n, function(i) {
    x[id == id.unq[i],, drop = FALSE][1,, drop = FALSE]}))
  x <- do.call(rbind, lapply(1:n, function(i) {x[id == id.unq[i],, drop = FALSE]}))
  if(any(is.na(X))){
    RemoveMiss.Obj <- RemoveMiss.Fun(X)
    X <- RemoveMiss.Obj$X
    if(!is.null(RemoveMiss.Obj$id.remove)){
      id.remove <- id.unq[RemoveMiss.Obj$id.remove]
      id.unq <- setdiff(id.unq,id.remove)
      n <- length(id.unq)
      tm <- unlist(lapply(1:n,function(i){ tm[id == id.unq[i] ]    }))
      y <- unlist(lapply(1:n,function(i){ y[id == id.unq[i] ]    }))
      x <- do.call(rbind, lapply(1:n, function(i) {x[id == id.unq[i],, drop = FALSE]}))
      id <- unlist(lapply(1:n,function(i){ id[id == id.unq[i] ]    }))
    }
  }
  p <- ncol(X)
  xvar.names <- colnames(X)
  if(family == "Continuous"){
    Ymean <- mean(y, na.rm = TRUE)
    Ysd <- sd(y, na.rm = TRUE)
    if (Ysd < 1e-6) {
      Ysd <- 1
    }
  }
  else{
    Ymean <- 0
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
  tm.list <- lapply(1:n, function(i) {tm[id == id.unq[i]]})
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
  bst.frac <- is.hidden.bst.frac(user.option)
  if(bootstrap == "by.user"){
    if(missing(bst.frac)){
      bst.frac <- 0.632
    }
    samp.mat <- matrix(NA,nrow = n,ncol = M)
    for(i in 1:M){
      samp.value <- (sample(1:n,floor(bst.frac*n),replace = FALSE))
      samp.value <- sort(c(samp.value,sample(samp.value, n - length(samp.value),replace = TRUE)))
      samp.value <- unlist(lapply(1:n,function(i){
        sum(samp.value == i)
      }))
      samp.mat[,i] <- samp.value
    }
  }
  if (ntree == 1) {
    nodesize <- max(1, round(n/(2 * K)))
    mtry <- df.D + p
  }
  else {
    nodedepth <- max(0, log(max(0, K), base = 2))
    nodesize <- 1
    mtry <- NULL
    if (missing(lambda) || lambda < 0) {
      lambda <- 0
    }
  }
  if (ntree > 1) {
    if (univariate) {
      learnerUsed <- "forest learner"
    }
    else {
      learnerUsed <- "mforest learner"
    }
  }
  else {
    if (df.D == 1) {
      learnerUsed <- "tree learner"
    }
    else {
      learnerUsed <- "mtree-Pspline learner"
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
  l_pred <- lapply(1:n, function(i) {rep(0, ni[i])})
  mu <- GetMu(Linear_Predictor = l_pred,Family = family)
  if (ntree == 1) {
    baselearner <- membership.list <- gamma.list <- vector("list", length = M)
  }
  else {
   membership.list <- gamma.list <- NULL
   baselearner <- vector("list", length = M)
  }
  if (!univariate) {
    lambda.vec <- phi.vec <- rho.vec <- rep(0, M)
  }
  lambda.initial <- Ysd^2
  rho.fit.flag <- TRUE
  rho.tree.grad <- 0
  rho.hide <- is.hidden.rho(user.option)
  if (!is.null(rho.hide) && (rho.hide >= 0 && rho.hide < 1)) {
    rho.fit.flag <- FALSE
    rho <- rho.hide
  }
  else {
    rho <- 0
  }
  sigma <- phi <- 1
  if (!lambda.est.flag) {
    sigma <- sigma.robust(lambda, rho)
  }
  Y.names <- paste("Y", 1:df.D, sep = "")
  rfsrc.f <- as.formula(paste("Multivar(", paste(Y.names, collapse = ","), paste(") ~ ."), sep = ""))
  cv.flag <- cv.flag && (ntree == 1)
  cv.lambda.flag <- cv.flag && is.hidden.CVlambda(user.option) && lambda.est.flag
  cv.rho.flag <- cv.flag && is.hidden.CVrho(user.option) && rho.fit.flag
  if (cv.flag) {
    mu.cv.list <- vector("list", M)
    l_pred.cv <- lapply(1:n, function(i) {rep(0, ni[i])})
    mu.cv <- GetMu(Linear_Predictor = l_pred.cv,Family = family)
    l_pred.i <- lapply(1:n, function(i) {
      lapply(1:n, function(j) {rep(0, ni[j])})
    })
    gamma.i.list <- lapply(1:M,function(m){
      vector("list", length = n)
    })
    err.rate <- matrix(NA, M, 2)
    colnames(err.rate) <- c("l1", "l2")
    if(family == "Continuous"){
      Ymean.i <- sapply(1:n, function(i) {
        mean(unlist(Yorg[-i]), na.rm = TRUE)
      })
      Ysd.i <- sapply(1:n, function(i) {
        sd.i <- sd(unlist(Yorg[-i]), na.rm = TRUE)
        if (sd.i < 1e-6) {
          1
        }
        else {
          sd.i
        }
      })
    }
    else{
      Ymean.i <- unlist(lapply(1:n,function(i){ Ymean  }))
      Ysd.i   <- unlist(lapply(1:n,function(i){ Ysd  }))
    }
  }
  else {
    err.rate <- rmse <- Mopt <- NULL
  }
  if(family == "Continuous"){
    NR.iter <- 1
  }
  if (verbose) pb <- txtProgressBar(min = 0, max = M, style = 3)
  for (m in 1:M) {
  if (verbose) setTxtProgressBar(pb, m)
    if (verbose && m == M) cat("\n")
    VMat <- lapply(1:n,function(i){
      VarTemp <- matrix(rho*phi,ni[i],ni[i])
      diag(VarTemp) <- phi
      VarTemp
    })
    inv.VMat <- lapply(1:n,function(i){
      out <- tryCatch({ qr.solve(VMat[[i]])},
                      error = function(ex){NULL})
      if(!is.null(out)){
        out <- out
      }else{
        out <- diag(phi,nrow(VMat[[i]]))
      }
      out
    })
    H_Mu <- Transform_H(Mu = mu, Family = family)
    if (mod.grad == FALSE) {
      gm.mod <- t(matrix(unlist(lapply(1:n, function(i) {
         t(D[[i]])%*%H_Mu[[i]]%*%inv.VMat[[i]]%*%(Y[[i]] - mu[[i]])  
      })), nrow = df.D))
    }
    else {
      gm.mod <- t(matrix(unlist(lapply(1:n, function(i) {
        t(D[[i]])%*%H_Mu[[i]]%*%(Y[[i]] - mu[[i]])
      })), nrow = df.D))
    }
    incoming.data <- cbind(gm.mod, X)
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
                         forest.wt = TRUE, 
                         memebership = TRUE)
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
                         bootstrap = bootstrap,
                         samp = if(bootstrap == "by.user") samp.mat[,m,drop = FALSE] else NULL,
                         membership = TRUE,
                         na.action = "na.impute",
                         nimpute = 1)
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
    }
    if (ntree == 1) {
      if (lambda.est.flag) {
        transf.data <- papply(1:n, function(i) {
          if (ni[i] > 1) {
            ci <- rho.inv.sqrt(ni[i], rho)
            R.inv.sqrt <- (diag(1, ni[i]) - matrix(ci, ni[i], ni[i])) / sqrt(1 - rho)
            V.inv.sqrt <- phi^(-1/2)*R.inv.sqrt
          }
          else {
            R.inv.sqrt <- cbind(1)
            V.inv.sqrt <- phi^(-1/2)*R.inv.sqrt
          }
          if (cv.lambda.flag) {
            Ynew <- V.inv.sqrt %*% (Y[[i]] - mu.cv[[i]])
          }
          else {
            Ynew <- V.inv.sqrt %*% (Y[[i]] - mu[[i]])
          }
            mu.2 <- GetMu_Lambda(Linear_Predictor = 2*l_pred[[i]],Family = family)
            LambdaD <- Transform_H(mu.2,Family = family)%*%D[[i]]
            Xnew <- V.inv.sqrt %*% LambdaD[, 1, drop = FALSE]
            Znew <- V.inv.sqrt %*% LambdaD[, -1, drop = FALSE] %*% pen.inv.sqrt.matx
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
        sigma <- sigma.robust(lambda, rho) 
      }
      gamma <- lapply(1:Kmax, function(k) {
        pt.k <- (membership == k)
        if (sum(pt.k) > 0) {
          which.pt.k <- which(pt.k == TRUE)
          seq.pt.k <- seq(length(which.pt.k) )
          gamma.NR.update <- rep(0,df.D)
          for(Iter in 1:NR.iter){
          mu.NR.update <- lapply(which.pt.k,function(i){
               l_pred_gamma <- l_pred[[i]]  + c(D[[i]]%*%gamma.NR.update)
               out <- GetMu(Linear_Predictor = l_pred_gamma, Family = family)
               out
          })
          CalD.i <- lapply(seq.pt.k,function(i){
            out_H_Mat <- Transform_H(Mu = mu.NR.update[[i]], Family = family)
            out <- out_H_Mat%*%D[[ which.pt.k[i] ]]
            out
          })
          HesMat.temp <- Reduce("+",lapply(seq.pt.k,function(i){
            t(CalD.i[[i]])%*%inv.VMat[[ which.pt.k[i] ]]%*%CalD.i[[i]]
          }))
          HesMat <- HesMat.temp + (lambda*pen.lsq.matx)
          ScoreVec.temp <-  Reduce("+",lapply(seq.pt.k,function(i){
            t(CalD.i[[i]])%*%inv.VMat[[ which.pt.k[i] ]]%*%(Y[[ which.pt.k[i]   ]] - mu.NR.update[[ i ]] )
          }))
          ScoreVec <- ScoreVec.temp - (lambda*(pen.lsq.matx%*%gamma.NR.update))
          qr.obj <- tryCatch({qr.solve(HesMat, ScoreVec)}, error = function(ex){NULL})
          if (!is.null(qr.obj)) {
            qr.obj <- qr.obj
          }
          else {
            qr.obj <- rep(0, df.D)
          }
          gamma.NR.update <- gamma.NR.update + qr.obj
          }
          gamma.NR.update
        }
        else {
          rep(0, df.D)
        }
      })
      gamma.matx <- matrix(0, Kmax, df.D + 1)
      gamma.matx[, 1] <- sort(unique(membership.org))
      gamma.matx[, 2:(df.D+1)] <- matrix(unlist(gamma), ncol = df.D, byrow = TRUE)
      gamma.list[[m]] <- gamma.matx
      bhat <- t(matrix(unlist(lapply(1:n, function(i) {
        gamma[[membership[i]]]})), nrow = df.D) * nu.vec)
      l_pred <- lapply(1:n, function(i) { l_pred[[i]] + c(D[[i]] %*% bhat[i, ])  })
      mu <- GetMu(Linear_Predictor = l_pred,Family = family)
      if (cv.flag) {
        oob <- which(rfsrc.obj$inbag == 0) 
        l_pred.i <- lapply(1:n,function(i) {
          if( any(i == oob)){
          mem.i <- membership[i]  
          l_pred.ij <- l_pred.i[[i]]
          gamma.i <- lapply(1:Kmax, function(k) {
            pt.k <- (membership == k)
            which.pt.k <- setdiff(which(pt.k == TRUE) ,i)
            if (sum(pt.k) > 0 && length(which.pt.k) > 0) {
              seq.pt.k <- seq(length(which.pt.k) )
              gamma.NR.update <- rep(0,df.D)
              for(Iter in 1:NR.iter){
                mu.NR.update <- lapply(which.pt.k,function(j){
                  l_pred_gamma <- l_pred.ij[[j]]  + c(D[[j]]%*%gamma.NR.update)
                  out <- GetMu(Linear_Predictor = l_pred_gamma, Family = family)
                  out
                })
                CalD.i <- lapply(seq.pt.k,function(j){
                  out_H_Mat <- Transform_H(Mu = mu.NR.update[[j]], Family = family)
                  out <- out_H_Mat%*%D[[ which.pt.k[j] ]]
                  out
                })
                HesMat.temp <- Reduce("+",lapply(seq.pt.k,function(j){
                  t(CalD.i[[j]])%*%inv.VMat[[ which.pt.k[j] ]]%*%CalD.i[[j]]
                }))
                HesMat <- HesMat.temp + (lambda*pen.lsq.matx)
                ScoreVec.temp <-  Reduce("+",lapply(seq.pt.k,function(j){
                  t(CalD.i[[j]])%*%inv.VMat[[ which.pt.k[j] ]]%*%(Y[[ which.pt.k[j]   ]] - mu.NR.update[[ j ]] )
                }))
                ScoreVec <- ScoreVec.temp - (lambda*(pen.lsq.matx%*%gamma.NR.update))
                qr.obj <- tryCatch({qr.solve(HesMat, ScoreVec)}, error = function(ex){NULL})
                if (!is.null(qr.obj)) {
                  qr.obj <- qr.obj
                }
                else {
                  qr.obj <- rep(0, df.D)
                }
                gamma.NR.update <- gamma.NR.update + qr.obj
              }
              gamma.NR.update
            }
            else {
              rep(0, df.D)
            }
          })
          gamma.matx.i <- matrix(0, Kmax, df.D + 1)
          gamma.matx.i[, 1] <- sort(unique(membership.org))
          gamma.matx.i[, 2:(df.D+1)] <- matrix(unlist(gamma.i), ncol = df.D, byrow = TRUE)
          gamma.i.list[[m]][[i]] <<- gamma.matx.i 
          lapply(1:n,function(j) {
            which.j <- which(gamma.matx.i[, 1] == membership.org[j])
            l_pred.i[[i]][[j]] + c(D[[j]] %*% (gamma.matx.i[which.j, -1] * nu.vec))
          })
      }
        else {
          l_pred.i[[i]]
      } 
      })   
        l_pred.cv <- lapply(1:n,function(i){l_pred.i[[i]][[i]]})
        mu.cv <- GetMu(Linear_Predictor = l_pred.cv,Family = family)
        l_pred.cv.org <- lapply(1:n,function(i){l_pred.cv[[i]] * Ysd.i[i] + Ymean.i[i]})
        mu.cv.org <- GetMu(Linear_Predictor = l_pred.cv.org,Family = family)
        mu.cv.list[[m]] <- mu.cv.org
        err.rate[m, ] <- c(l1Dist(Yorg, mu.cv.org), l2Dist(Yorg, mu.cv.org))
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
          YnewSum <- colSums(fwt.i[pt.i] * gm.mod[pt.i,, drop = FALSE])
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
      bhat <- t(bhat * nu.vec)
      mu <- lapply(1:n, function(i) {mu[[i]] + D[[i]] %*% bhat[i, ]})
      baselearner[[m]] <- c(baselearner[[m]],
                            list(gm = gm.mod),
                            list(Xnew = Xnew),
                            list(pen = sigma * pen.lsq.matx))
    }
    if (!univariate && rho.fit.flag) {
      if (cv.rho.flag) {
        resid.data <- data.frame(y  = unlist(lapply(1:n, function(i) {Y[[i]] - mu.cv[[i]]})),
                                 x,
                                 tm = unlist(lapply(1:n, function(i) {tm[id == id.unq[i]]})),
                                 id = unlist(lapply(1:n, function(i) {rep(id.unq[i], ni[i])})))
      }
      else {
        resid.data <- data.frame(y  = unlist(lapply(1:n, function(i) {Y[[i]] - mu[[i]]})),
                                 x,
                                 tm = unlist(lapply(1:n, function(i) {tm[id == id.unq[i]]})),
                                 id = unlist(lapply(1:n, function(i) {rep(id.unq[i], ni[i])})))
      }
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
        rho <- as.numeric(coef(gls.obj$modelStruct$corStruc, unconstrained = FALSE))
        rho <- max(min(0.999, rho, na.rm = TRUE), -0.999)
      }
    }
    if (!univariate) {
      phi.vec[m] <- phi * Ysd^2
      rho.vec[m] <- rho
    }
    if (!univariate) {
      sigma <- sigma.robust(lambda, rho)
      lambda.vec[m] <- lambda
    }
  }
  l_pred <- lapply(1:n, function(i) {c(l_pred[[i]] * Ysd + Ymean)})
  mu <- GetMu(Linear_Predictor = l_pred,Family = family)
  y <- lapply(1:n, function(i) {y[id == id.unq[i]]})
  if (cv.flag) {
    diff.err <- abs(err.rate[, "l2"] - min(err.rate[, "l2"], na.rm = TRUE))
    diff.err[is.na(diff.err)] <- 1
      if (sum(diff.err < Ysd * eps) > 0) {
        Mopt <- min(which(diff.err < eps))
      }
      else {
        Mopt <- M
      }
    rmse <- err.rate[Mopt, "l2"]
    mu <- lapply(1:n,function(i){mu.cv.list[[Mopt]][[i]]})
  }
  obj <- list(x = X,
              xvar.names = xvar.names,
              time = lapply(1:n, function(i) {tm[id == id.unq[i]]}),
              id = id,
              y = y,
              family = family,
              ymean = Ymean,
              ysd = Ysd,
              n = n,
              ni = ni,
              tm.unq = tm.unq,
              gamma = gamma.list,
              mu = mu,
              lambda = lambda.vec,
              phi = phi.vec,
              rho = rho.vec,
              baselearner = baselearner,
              membership = membership.list,
              X.tm = X.tm,
              D = D,
              d = d,
              pen.ord = pen.ord,
              K = K,
              M = M,
              nu = nu,
              ntree = ntree,
              cv.flag = cv.flag,
              err.rate = if (!is.null(err.rate)) err.rate / Ysd else NULL,
              rmse = if (!is.null(rmse)) rmse / Ysd else NULL,
              Mopt = Mopt,
              gamma.i.list = if(cv.flag) gamma.i.list else NULL,
              forest.tol = forest.tol)
  class(obj) <- c("boostmtree", "grow", learnerUsed)
  invisible(obj)
}
