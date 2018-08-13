####**********************************************************************
####**********************************************************************
####
####  BOOSTED MULTIVARIATE TREES FOR LONGITUDINAL DATA (BOOSTMTREE)
####  Version 1.3.0 (_PROJECT_BUILD_ID_)
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


marginalPlot <- function (obj,
                         xvar.names,
                         tm,
                         subset,
                         plot.it = TRUE,
                         ...)
{
  if (sum(inherits(obj, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }
  if (missing(xvar.names)) {
    xvar.names <- colnames(obj$x)
  }
  xvar.names <- intersect(xvar.names, colnames(obj$x))
  if (length(xvar.names) == 0) {
    stop("x-variable names provided do not match original variable names")
  }
  n.xvar <- length(xvar.names)
  tmOrg <- sort(unique(unlist(obj$time)))
  if (missing(tm)) {
    tm.q <- unique(quantile(tmOrg, (1:9)/10, na.rm = TRUE))
    tm.pt <- sapply(tm.q, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  else {
    tm.q <- tm
    tm.pt <- sapply(tm, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  n.tm <- length(tm.pt)
  if (!missing(subset)) {
    obj$x <- obj$x[subset,, drop = FALSE]
  }
  n <- nrow(obj$x)
  if(n.tm == 1){
      muhat <- cbind(matrix(unlist(predict.boostmtree(object = obj,importance = FALSE)$muhat),nrow = n,byrow = TRUE)[,tm.pt])
    }else
      {
      muhat <- matrix(unlist(predict.boostmtree(object = obj,importance = FALSE)$muhat),nrow = n,byrow = TRUE)[,tm.pt]
    }
  lo.obj <- lapply(1:n.xvar, function(nm){
    x <- obj$x[, xvar.names[nm]]
    lo.fit <- lapply(1:n.tm,function(nt){
      fit <- lowess(x,muhat[,nt])
      cbind(fit$x,fit$y)
    })
    names(lo.fit) <- paste("time = ",tm.q,sep="")
    lo.fit
  })
  names(lo.obj) <- xvar.names
  if(plot.it){
    if(n.xvar > 1){
      pdf(file = "MarginalPlot.pdf",width = 10,height = 10)
    }
    for(pp in 1:n.xvar){
      xmin <- min(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,1]   })))
      xmax <- max(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,1]   })))
      ymin <- min(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,2]   })))
      ymax <- max(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,2]   })))
      plot(lo.obj[[pp]][[1]][,1],lo.obj[[pp]][[1]][,2],type = "n",xlim=c(xmin,xmax),ylim=c(ymin,ymax) ,
           xlab = xvar.names[pp],ylab = "Predicted response")
      for(nn in 1:n.tm){
      lines(lo.obj[[pp]][[nn]][,1],lo.obj[[pp]][[nn]][,2],type = "l",col = nn)
      }
    }
    if(n.xvar > 1){
    dev.off()
    print(paste("Plot is stored in the directory:",getwd(),sep=" "))
    }
  }
  return(invisible(lo.obj))
}
