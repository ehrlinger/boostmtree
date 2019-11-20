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


vimpPlot <- function(vimp,
                     Time_Interaction = TRUE,
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
                     seplim = NULL,
                     eps = 0.1,
                     Width_Bar = 1
)
{
  if(is.null(vimp) ){
    stop("vimp is not present in the object")
  }
  if(Time_Interaction){
    p <- floor(length(vimp)/2)
  }else
  {
    p <- length(vimp)  
  }
  if(is.null(xvar.names)){
    xvar.names <- paste("x",1:p,sep="")
  }
  vimp <- vimp*100
  if(!Time_Interaction){
    ylim <- range(vimp) + c(-0,ymaxlim)
    yaxs <- pretty(ylim)
    yat <- abs(yaxs)
    bp <- barplot(pmax(as.matrix(vimp),0),beside=T,width = Width_Bar,col=col,ylim=ylim,yaxt="n",main = main,cex.lab=cex.lab)
    text(c(bp), pmax(as.matrix(vimp),0) + eps, rep(xvar.names, 3),srt=90,adj= 0.0,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
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
      xaxishead <- c(floor(p/4),floor(p/4))
    }
    bp1 <- barplot(pmax(as.matrix(vimp.x),0),width = Width_Bar,horiz = FALSE,beside=T,col=col,ylim=ylim,yaxt="n",ylab = ylbl,cex.lab=cex.lab,
                   main = main,line = 2)
    text(c(bp1), pmax(as.matrix(vimp.x),0) + eps, rep(xvar.names, 3),srt=90,adj= 0.0,cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[2],yaxishead[2],labels = subhead.labels[2],cex = subhead.cexval)
    bp2 <- barplot(-pmax(as.matrix(vimp.time),0) - eps,width = Width_Bar,horiz = FALSE,beside=T,col=col,add=TRUE,yaxt="n")
    #text(c(bp2), -4, rep(xvar.names, 3),srt=270,adj= 0,yaxt="n",cex=if(!is.null(cex.xlab)) cex.xlab else 1)
    text(xaxishead[1],-yaxishead[1],labels = subhead.labels[1],cex = subhead.cexval)
    axis(2,yaxs,yat)
  }
}
