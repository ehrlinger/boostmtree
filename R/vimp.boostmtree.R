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


vimp.boostmtree <- function(object,x.names = NULL,joint = FALSE){
  if (sum(inherits(object, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(boostmtree, grow)' or '(boostmtree, predict)'")
  }
  if (sum(inherits(object, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {
  if(!object$cv.flag){
    stop("The grow object of boostmtree does not include in-sample CV estimates")
  }  
  X <- object$x
  P <- ncol(X)
  xvar.names <- object$xvar.names
  if(is.null(x.names)){
    vimp_set <- 1:P
    x_Names <- xvar.names
  }else
  {
      n.x.names <- length(x.names)
      vimp_set <- match(x.names,xvar.names)
      if(any(is.na( vimp_set ))){
        stop("x.names do not match with variable names from original data")
      }
      if(joint){
        x_Names <- "joint_vimp"
      }else
      {
        x_Names <- x.names
      }
  }
  if(joint){
    p <- 1
  }else
  {
    p <- length(vimp_set)  
  }
  Y <- object$y
  Ymean <- object$ymean
  Ysd <- object$ysd
  df.D <- ncol(object$X.tm)
  Mopt <- object$Mopt
  ni <- object$ni
  K <- object$K
  gamma.i.list <- object$gamma.i.list
  membership <- object$membership
  D <- object$D
  n <- object$n
  rmse <- object$rmse*object$ysd
  nu.vec <- rep(object$nu[1],df.D)
  family <- object$family
  oob.list <- vector("list",Mopt)
  membershipNoise.list <- lapply(1:Mopt,function(m){
    oob <- which(object$baselearner[[m]]$inbag == 0)
    oob.list[[m]] <<- oob
    n.oob <- length(oob)
    Xnoise <- do.call(rbind, lapply(1:p, function(k) {
      X.k <- X[oob,, drop = FALSE]
      if(joint){
        X.k[, vimp_set] <- X.k[sample(nrow(X.k)), vimp_set]
      }else
      {
        X.k[, vimp_set[k] ] <- sample(X.k[, vimp_set[k] ])
      }
      X.k
    }))
    membershipNoise <- c(predict.rfsrc(object$baselearner[[m]],
                                       newdata = Xnoise,
                                       membership = TRUE,
                                       ptn.count = K,
                                       importance = "none",
                                       na.action = "na.impute")$ptn.membership)
    membershipNoise <- matrix(membershipNoise,nrow = n.oob,byrow = FALSE)
    membershipNoise
  })
  if(df.D > 1){
    vimp.main <- vimp.int <-rep(NA,p)
    for(k in 1:p){
      l_pred.vimp <- lapply(1:n,function(i){
        l_pred.main.i <- l_pred.int.i <- rep(0,ni[i])
        if(k == p){ l_pred.time <- rep(0,ni[i])  }
        NullObj <- lapply(1:Mopt,function(m){
          if( any(i == oob.list[[m]] )){
            membershipNoise.i <- membershipNoise.list[[m]][ which(oob.list[[m]] == i) , k ,drop = TRUE]  
            membershipOrg.i.vec <- gamma.i.list[[m]][[i]][,1,drop = TRUE]
            gamma.noise.i <- t(gamma.i.list[[m]][[i]][which(membershipOrg.i.vec == membershipNoise.i),-1,drop = FALSE])
            membershipOrg.i <- membership[[m]][ i ]
            gamma.org.i <- t(gamma.i.list[[m]][[i]][which(membershipOrg.i.vec == membershipOrg.i),-1,drop = FALSE])
            gamma.main <- cbind(c(gamma.noise.i[1,1],gamma.org.i[-1,1]))
            out.main <- c(D[[i]] %*% (gamma.main * nu.vec))
            gamma.int <- cbind(c(gamma.org.i[1,1],gamma.noise.i[-1,1]))
            out.int <- c(D[[i]] %*% (gamma.int * nu.vec))
            if(k == p){
              n.D <- nrow(D[[i]])
              out.time <- c(D[[i]][sample(1:n.D,n.D,replace = TRUE),,drop = FALSE] %*% (gamma.org.i * nu.vec))
            }
          }
          else{
            out.main <- out.int <- rep(0,ni[i])
            if(k == p){ out.time <- rep(0,ni[i])  }
          }
          l_pred.main.i <<- l_pred.main.i + out.main
          l_pred.int.i <<- l_pred.int.i + out.int
          if(k == p){ l_pred.time <<- l_pred.time + out.time }
          NULL
        })
        list(l_pred.main = l_pred.main.i,l_pred.int = l_pred.int.i,l_pred.time = if (k == p) l_pred.time else NULL  )
      })
      mu.main <- lapply(1:n,function(i){  GetMu(Linear_Predictor = l_pred.vimp[[i]]$l_pred.main * Ysd + Ymean ,Family = family)   })
      mu.int <- lapply(1:n,function(i){  GetMu(Linear_Predictor = l_pred.vimp[[i]]$l_pred.int * Ysd + Ymean,Family = family)   })
      err.rate.main <- l2Dist(Y, mu.main)
      vimp.main[k] <-  err.rate.main - rmse   
      err.rate.int <- l2Dist(Y, mu.int)
      vimp.int[k] <-  err.rate.int - rmse   
      if(k == p){  
        mu.time <- lapply(1:n,function(i){ GetMu(Linear_Predictor = l_pred.vimp[[i]]$l_pred.time * Ysd + Ymean,Family = family)   })
        err.rate.time <- l2Dist(Y, mu.time)
        vimp.time <-  err.rate.time - rmse  
      }
    }
    names(vimp.main) <- x_Names
    names(vimp.int) <- paste(x_Names, "time", sep=":")
    names(vimp.time) <- "time"
    vimp <- c(vimp.main,vimp.int,vimp.time)  
  }
  else{
    vimp <- rep(NA,p)
    for(k in 1:p){
      l_pred.k <- lapply(1:n,function(i){
        l_pred.i <- rep(0,ni[i])    
        NullObj <- lapply(1:Mopt,function(m){
          if( any(i == oob.list[[m]] )){
            membershipNoise.i <- membershipNoise.list[[m]][ which(oob.list[[m]] == i)   , k ,drop = TRUE]  
            membershipOrg.i.vec <- gamma.i.list[[m]][[i]][,1,drop = TRUE]
            gamma.noise.i <- t(gamma.i.list[[m]][[i]][which(membershipOrg.i.vec == membershipNoise.i),-1,drop = FALSE])
            out <- c(D[[i]] %*% (gamma.noise.i * nu.vec))
          }
          else
          {
            out <- rep(0,ni[i])
          }
          l_pred.i <<- l_pred.i + out
          NULL
        })
        l_pred.i
      })
      mu.k <- lapply(1:n,function(i){ GetMu(Linear_Predictor = l_pred.k[[i]] * Ysd + Ymean , Family = family) })
      err.rate.k <- l2Dist(Y, mu.k)
      vimp[k] <- err.rate.k - rmse 
    }
    names(vimp) <- x_Names
   }
}
  else {
    X <- object$x
    P <- ncol(X)
    xvar.names <- object$xvar.names
    if(is.null(x.names)){
      vimp_set <- 1:P
      x_Names <- xvar.names
    }else
    {
      n.x.names <- length(x.names)
      vimp_set <- match(x.names,xvar.names)
      if(any(is.na( vimp_set ))){
        stop("x.names do not match with variable names from original data")
      }
      if(joint){
        x_Names <- "joint_vimp"
      }else
      {
        x_Names <- x.names
      }
    }
    if(joint){
      p <- 1
    }else
    {
      p <- length(vimp_set)  
    }
    Y <- object$y
    if(is.null(Y)){
      stop("Response is not provied in the predict object")
    }
    Ymean <- object$ymean
    Ysd <- object$ysd
    n <- object$n
    K <- object$K
    ni <- object$ni
    df.D <- object$df.D
    D <- object$D
    nu.vec <- object$nu.vec
    Mopt  <- object$Mopt
    baselearner <- object$baselearner
    gamma <- object$gamma
    membership <- object$membership
    rmse <- object$rmse*Ysd
    family <- object$family
    membershipNoise <- lapply(1:Mopt, function(m) {
      Xnoise <- do.call(rbind, lapply(1:p, function(k) {
        X.k <- X
        if(joint){
          X.k[, vimp_set] <- X.k[sample(nrow(X.k)), vimp_set]
        }else
        {
          X.k[, vimp_set[k] ] <- sample(X.k[, vimp_set[k] ])
        }
        X.k
      }))
      options(rf.cores = 1, mc.cores = 1)
      c(predict.rfsrc(baselearner[[m]],
                      newdata = Xnoise,
                      membership = TRUE,
                      ptn.count = K,
                      na.action = "na.impute",
                      importance = "none")$ptn.membership)
    })
    if(df.D > 1){
      vimp_main <- vimp_int <- rep(NA,p)
      vimp_time <- NA
      nullObj <- lapply(1:p, function(k) {
        l_pred_vimp <- lapply(1:n,function(i){
          l_pred_main.i <- l_pred_int.i <- l_pred_time.i <- rep(0,ni[i])
          NullObj <- lapply(1:Mopt,function(m){
            orgMembership  <- gamma[[m]][, 1]
            gamma.Org      <- gamma[[m]][match(membership[[m]][i], orgMembership), -1, drop = FALSE]
            membership.k   <- membershipNoise[[m]][((k-1) * n + 1):(k * n)]
            membership.k.i <- membership.k[i]
            gamma.Noise    <- gamma[[m]][match(membership.k.i, orgMembership), -1, drop = FALSE]
            gamma.main     <- cbind(c(gamma.Noise[1],gamma.Org[-1]))
            gamma.int      <- cbind(c(gamma.Org[1],gamma.Noise[-1]))
            l_pred_main.i  <<- l_pred_main.i + c(D[[i]]%*%(gamma.main*nu.vec))
            l_pred_int.i   <<- l_pred_int.i  + c(D[[i]]%*%(gamma.int*nu.vec))
            if(k == p){
              n.D <- nrow(D[[i]])
              l_pred_time.i <<- l_pred_time.i + D[[i]][sample(1:n.D,n.D,replace = TRUE),,drop = FALSE]%*%t(gamma.Org*nu.vec)
            }
            NULL
          })
          list(l_pred_main.i,l_pred_int.i,l_pred_time.i)
        })
        mu_main <- lapply(1:n,function(i){ GetMu(Linear_Predictor = l_pred_vimp[[i]][[1]] * Ysd + Ymean, Family = family) })
        mu_int <- lapply(1:n,function(i){ GetMu(Linear_Predictor = l_pred_vimp[[i]][[2]] * Ysd + Ymean, Family = family) })
        err.rate.main  <- l2Dist(Y, mu_main)
        err.rate.int   <- l2Dist(Y, mu_int)
        vimp_main[k]  <<- err.rate.main - rmse
        vimp_int[k]   <<- err.rate.int - rmse
        if(k == p){
          mu_time <- lapply(1:n,function(i){ GetMu(Linear_Predictor = l_pred_vimp[[i]][[3]] * Ysd + Ymean, Family = family) })
          err.rate.time  <- l2Dist(Y, mu_time)
          vimp_time     <<- err.rate.time - rmse
        }
        NULL
      })
      names(vimp_main) <- x_Names
      names(vimp_int) <- paste(x_Names, "time", sep=":")
      names(vimp_time) <- "time"
      vimp <- c(vimp_main,vimp_int,vimp_time)  
    }
    else {
      vimp <- rep(NA,p)
      nullObj <- lapply(1:p, function(k) {
        l_pred_vimp <- lapply(1:n,function(i){
          l_pred.i <- rep(0,ni[i])
          NullObj <- lapply(1:Mopt,function(m){
            orgMembership  <- gamma[[m]][, 1]
            membership.k   <- membershipNoise[[m]][((k-1) * n + 1):(k * n)]
            membership.k.i <- membership.k[i]
            gamma.Noise    <- gamma[[m]][match(membership.k.i, orgMembership), -1, drop = FALSE]
            l_pred.i       <<- l_pred.i  + c(D[[i]]%*%t(gamma.Noise*nu.vec))
            NULL
          })
          l_pred.i
        })
        mu_vimp <- lapply(1:n,function(i){ GetMu(Linear_Predictor = l_pred_vimp[[i]] * Ysd + Ymean, Family = family)  })
        err.rate  <- l2Dist(Y, mu_vimp)
        vimp[k]       <<- err.rate - rmse
        NULL
      })
      names(vimp) <- xvar.names
    }
  }
  obj <- vimp/rmse
  invisible(obj)
}
