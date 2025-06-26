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










#' Variable Importance (VIMP) plot
#' 
#' Barplot displaying VIMP.
#' 
#' Barplot displaying VIMP. If the analysis is for the univariate case, VIMP is
#' displayed above the x-axis. If the analysis is for the longitudinal case,
#' VIMP for covariates (main effects) are shown above the x-axis while VIMP for
#' covariate-time interactions (time interaction effects) are shown below the
#' x-axis. In either case, negative vimp value is set to zero.
#' 
#' @param vimp VIMP values.
#' @param Q_set Provide names for various levels of nominal or ordinal
#' response.
#' @param Time_Interaction Whether VIMP is estimated from a longitudinal data,
#' in which case VIMP is available for covariate and covariate-time
#' interaction. Default is TRUE. If FALSE, VIMP is assumed to be estimated from
#' a cross-sectional data.
#' @param xvar.names Names of the covariates. If NULL, names are assigned as
#' x1, x2,...,xp.
#' @param cex.xlab Magnification of the names of the covariates above (and
#' below) the barplot.
#' @param ymaxlim By default, we use the range of the vimp values for the
#' covariates for the ylim. If one wants to extend the ylim, add the amount
#' with which the ylim will extend above.
#' @param ymaxtimelim By default, we use the range of the vimp values for the
#' covariates-time for the ylim. If one wants to extend the ylim, add the
#' amount with which the ylim will extend below. Argument only works for the
#' longitudinal setting.
#' @param subhead.cexval Magnification of the \code{subhead.labels}. Argument
#' only works for the longitudinal setting.
#' @param yaxishead This represent a vector with two values which are points on
#' the y-axis. Corresponding to the values, the lables for
#' \code{subhead.labels} is shown. First argument corresponds to covariate-time
#' interaction, whereas second argument is for the main effect. Argument only
#' works for the longitudinal setting.
#' @param xaxishead This represent a vector with two values which are points on
#' the x-axis. Corresponding to the values, the lables for
#' \code{subhead.labels} is shown. First argument corresponds to covariate-time
#' interaction, whereas second argument is for the main effect. Argument only
#' works for the longitudinal setting.
#' @param main Main title for the plot.
#' @param col Color of the plot.
#' @param cex.lab Magnification of the x and y lables.
#' @param subhead.labels Labels corresponding to the plot. Default is
#' "Time-Interactions Effects" for the barplot below x-axis, and "Main Effects"
#' for the barplot above x-axis.
#' @param ylbl Should labels for the sub-headings be shown on left side of the
#' y-axis.
#' @param seplim if \code{ylbl} is \code{TRUE}, the distance between the lables
#' of the sub-headings.
#' @param eps Amount of gap between the top of the barplot and variable names.
#' @param Width_Bar Width of the barplot.
#' @param path_saveplot Provide the location where plot should be saved. By
#' default the plot will be saved at temporary folder.
#' @param Verbose Display the path where the plot is saved?
#' @author Hemant Ishwaran, Amol Pande and Udaya B. Kogalur
#' @keywords plot
#' @examples
#' 
#' \dontrun{
#' ##------------------------------------------------------------
#' ## Synthetic example
#' ## high correlation, quadratic time with quadratic interaction
#' ##-------------------------------------------------------------
#' #simulate the data
#' dta <- simLong(n = 50, N = 5, rho =.80, model = 2,family = "Continuous")$dtaL
#' 
#' #basic boosting call
#' boost.grow <- boostmtree(dta$features, dta$time, dta$id, dta$y, 
#'               family = "Continuous",M = 300, cv.flag = TRUE)
#' vimp.grow <- vimp.boostmtree(object = boost.grow)              
#' 
#' # VIMP plot
#' vimpPlot(vimp = vimp.grow, ymaxlim = 20, ymaxtimelim = 20,
#'          xaxishead = c(3,3), yaxishead = c(65,65),
#'          cex.xlab = 1, subhead.cexval = 1.2)
#' }
#' 
vimpPlot <- function(vimp,
                     Q_set = NULL,
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
                     subhead.labels = c("Time-Interactions Effects", "Main Effects"),
                     ylbl = FALSE,
                     seplim = NULL,
                     eps = 0.1,
                     Width_Bar = 1,
                     path_saveplot = NULL,
                     Verbose = TRUE) {
  if (is.null(vimp)) {
    stop("vimp is not present in the object")
  }
  
  
  p <- nrow(vimp[[1]])
  if (is.null(p)) {
    p <- nrow(vimp)
  }
  if (is.null(xvar.names)) {
    xvar.names <- paste("x", 1:p, sep = "")
  }
  vimp <- lapply(vimp, function(v) {
    v * 100
  })
  n.Q <- ncol(vimp[[1]])
  if (is.null(Q_set)) {
    Q_set <- paste("V", seq(n.Q), sep = "")
  }
  for (q in 1:n.Q) {
    if (is.null(path_saveplot)) {
      path_saveplot <- tempdir()
    }
    Plot_Name <- if (n.Q == 1)
      "VIMPplot.pdf"
    else
      paste("VIMPplot_Prob(y = ", Q_set[q], ")", ".pdf", sep = "")
    pdf(
      file = paste(path_saveplot, "/", Plot_Name, sep = ""),
      width = 10,
      height = 10
    )
    if (!Time_Interaction) {
      ylim <- range(vimp[[1]][, q]) + c(-0, ymaxlim)
      yaxs <- pretty(ylim)
      yat <- abs(yaxs)
      bp <- barplot(
        pmax(as.matrix(vimp[[1]][, q]), 0),
        beside = TRUE,
        width = Width_Bar,
        col = col,
        ylim = ylim,
        yaxt = "n",
        main = main,
        cex.lab = cex.lab
      )
      text(
        c(bp),
        pmax(as.matrix(vimp[[1]][, q]), 0) + eps,
        rep(xvar.names, 3),
        srt = 90,
        adj = 0.0,
        cex = if (!is.null(cex.xlab))
          cex.xlab
        else
          1
      )
      axis(2, yaxs, yat)
    } else {
      vimp.x <- vimp[[1]][, q]
      vimp.time <- vimp[[2]][, q]
      ylim <- max(c(vimp.x, vimp.time)) * c(-1, 1) + c(-ymaxtimelim, ymaxlim)
      if (ylbl) {
        ylabel <- paste("Time-Interactions",
                        "Main Effects",
                        sep = if (!is.null(seplim))
                          seplim
                        else
                          "                   ")
      } else {
        ylabel <- ""
      }
      yaxs <- pretty(ylim)
      yat <- abs(yaxs)
      if (is.null(yaxishead)) {
        yaxishead <- c(-ylim[1], ylim[2])
      }
      if (is.null(xaxishead)) {
        xaxishead <- c(floor(p / 4), floor(p / 4))
      }
      bp1 <- barplot(
        pmax(as.matrix(vimp.x), 0),
        width = Width_Bar,
        horiz = FALSE,
        beside = TRUE,
        col = col,
        ylim = ylim,
        yaxt = "n",
        ylab = ylabel,
        cex.lab = cex.lab,
        main = main
      )
      text(
        c(bp1),
        pmax(as.matrix(vimp.x), 0) + eps,
        rep(xvar.names, 3),
        srt = 90,
        adj = 0.0,
        cex = if (!is.null(cex.xlab))
          cex.xlab
        else
          1
      )
      text(xaxishead[2],
           yaxishead[2],
           labels = subhead.labels[2],
           cex = subhead.cexval)
      bp2 <- barplot(
        -pmax(as.matrix(vimp.time), 0) - eps,
        width = Width_Bar,
        horiz = FALSE,
        beside = TRUE,
        col = col,
        add = TRUE,
        yaxt = "n"
      )
      # text(c(bp2), -4, rep(xvar.names, 3),srt=270,adj= 0,
      # yaxt="n",cex=if(!is.null(cex.xlab)) cex.xlab else 1)
      text(xaxishead[1],-yaxishead[1],
           labels = subhead.labels[1],
           cex = subhead.cexval)
      axis(2, yaxs, yat)
    }
    dev.off()
    if (Verbose) {
      cat("Plot will be saved at:", path_saveplot, sep = "", "\n")
    }
  }
}
