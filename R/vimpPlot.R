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


vimpPlot <- function(object,
                     xvar.names = NULL,
                     cex.xlab = NULL,
                     ymaxlim = 0,
                     ymaxtimelim = 0,
                     subhead.cexval = 1,
                     yaxishead = NULL,
                     xaxishead = NULL,
                     main = "Variable Importance (%)",
                     col = grey(.80),
                     cex.lab = 1.5,
                     subhead.labels = c("Time-Interactions Effects","Main Effects"),
                     ylbl = FALSE,
                     seplim = NULL
)
{
  if(is.null(object$vimp) ){
    stop("vimp is not present in the object")
  }
  vimp <- object$vimp
  if(is.null(xvar.names)){
    xvar.names <- colnames(object$x)
  }
  p <- ncol(object$x)
  n.vimp <- length(vimp)
  if(n.vimp == p ){
    univariate <- TRUE
  }else
  {
    univariate <- FALSE
  }
  if(univariate){
    vimp <- vimp*100
  }else
  {
    vimp <- (vimp[-n.vimp])*100
    n.vimp <- length(vimp)
  }
  if(univariate){
    ylim <- range(vimp) + c(0,ymaxlim)
    yaxs <- pretty(ylim)
    yat <- abs(yaxs)
    bp <- barplot(as.matrix(vimp),beside=T,col=col,ylim=ylim,yaxt="n",main = main,cex.lab=cex.lab)
    text(c(bp), pmax(as.matrix(vimp),0), rep(xvar.names, 3),srt=90,adj=-0.5,cex= if(!is.null(cex.xlab)) cex.xlab else 1 )
    axis(2,yaxs,yat)
  }else
  {
    vimp.x <- vimp[1:p]
    vimp.time <- vimp[-c(1:p)]
    ylim <- max(c(vimp.x,vimp.time)) * c(-1, 1) + c(-ymaxtimelim,ymaxlim)
    if(ylbl){
      ylbl <- paste("Time-Interactions", "Main Effects", sep = if(!is.null(seplim)) seplim else "                   " )
    }else
    {
      ylbl <- NULL
    }
    yaxs <- pretty(ylim)
    yat <- abs(yaxs)
    if(is.null(yaxishead)){
      yaxishead <- c(-ylim[1],ylim[2])
    }
    if(is.null(xaxishead)){
      xaxishead <- c(floor(n.vimp/4),floor(n.vimp/4))
    }
    bp1 <- barplot(pmax(as.matrix(vimp.x),0),beside=T,col=col,ylim=ylim,yaxt="n",ylab = ylbl,cex.lab=cex.lab,
                   main = main)
    text(c(bp1), pmax(as.matrix(vimp.x),0), rep(xvar.names, 3),srt=90,adj=-0.5,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[2],yaxishead[2],labels = subhead.labels[2],cex = subhead.cexval)
    bp2 <- barplot(-pmax(as.matrix(vimp.time),0),beside=T,col=col,add=TRUE,yaxt="n")
    text(c(bp2), -pmax(as.matrix(vimp.time),0), rep(xvar.names, 3),srt=90,adj=1.5,yaxt="n",cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[1],-yaxishead[1],labels = subhead.labels[1],cex = subhead.cexval)
    axis(2,yaxs,yat)
  }
}
