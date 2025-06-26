
#' Atrial Fibrillation Data
#' 
#' Atrial Fibrillation (AF) data is obtained from a randomized trial to study
#' the effect of surgical ablation as a treatment option for AF among patients
#' with persistent and long-standing persistent AF who requires mitral valve
#' surgery. Patients were randomized into two groups: mitral valve surgery with
#' ablation and mitral valve surgery without ablation. Patients in the ablation
#' group were further randomized into two types of procedure: pulmonary vain
#' isolation (PVI) and biatrial maze procedure. These patients were followed
#' weekly for a period of 12 months. The primary outcome of the study is the
#' presence/absence of AF (binary longitudinal response). Data includes 228
#' patients. From 228 patients, 7949 AF measurements are available with average
#' of 35 measurements per patient.
#' 
#' 
#' @name AF
#' @docType data
#' @format A list containing four elements: \enumerate{ \item The 84 patient
#' variables (features).  \item Time points (time).  \item Unique patient
#' identifier (id).  \item Presence or absence of AF (y).  }
#' @references Gillinov A. M., Gelijns A.C., Parides M.K., DeRose J.J.Jr.,
#' Moskowitz~A.J. et al. Surgical ablation of atrial fibrillation during mitral
#' valve surgery. \emph{The New England Journal of Medicine}
#' 372(15):1399--1408, 2015.
#' @keywords datasets
#' @examples
#' data(AF, package = "boostmtree")
NULL





#' Boosted multivariate trees for longitudinal data.
#' 
#' Multivariate extension of Friedman's (2001) gradient descent boosting method
#' for modeling longitudinal response using multivariate tree base learners.
#' Longitudinal response could be continuous, binary, nominal or ordinal.
#' Covariate-time interactions are modeled using penalized B-splines
#' (P-splines) with estimated adaptive smoothing parameter.
#' 
#' 
#' @name boostmtree-package
#' @docType package
#' @section Package Overview: This package contains many useful functions and
#' users should read the help file in its entirety for details.  However, we
#' briefly mention several key functions that may make it easier to navigate
#' and understand the layout of the package.
#' 
#' \enumerate{ \item \command{\link{boostmtree}}
#' 
#' This is the main entry point to the package.  It grows a multivariate tree
#' using user supplied training data.  Trees are grown using the
#' \pkg{randomForestSRC} R-package.
#' 
#' \item \command{\link{predict.boostmtree}} (\command{predict})
#' 
#' Used for prediction.  Predicted values are obtained by dropping the user
#' supplied test data down the grow forest.  The resulting object has class
#' \code{(rfsrc, predict)}.
#' 
#' }
#' @author Hemant Ishwaran, Amol Pande and Udaya B. Kogalur
#' @seealso \command{\link{partialPlot}}, \command{\link{plot.boostmtree}},
#' \command{\link{predict.boostmtree}}, \command{\link{print.boostmtree}},
#' \command{\link{simLong}}
#' @references Friedman J.H. (2001). Greedy function approximation: a gradient
#' boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.
#' 
#' Friedman J.H. (2002). Stochastic gradient boosting.  \emph{Comp. Statist.
#' Data Anal.}, 38(4):367--378.
#' 
#' Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B., Blackstone E.H.,
#' Ishwaran H. (2017).  Boosted multivariate trees for longitudinal data,
#' \emph{Machine Learning}, 106(2): 277--305.
#' @keywords package
NULL





#' Spirometry Data
#' 
#' Data consists of 9471 longitudinal evaluations of forced 1-second expiratory
#' volume (FEV1-percentage of predicted) after lung transplant from 509
#' patients who underwent lung transplant (LTx) at the Cleveland Clinic.
#' Twenty three patient/procedure variables were collected at the time of the
#' transplant.  The major objectives are to evaluate the temporal trend of FEV1
#' after LTx, and to identify factors associated with post-LTx FEV1 and
#' assessing the differences in the trends after Single LTx versus Double LTx.
#' 
#' 
#' @name spirometry
#' @docType data
#' @format A list containing four elements: \enumerate{ \item The 23 patient
#' variables (features).  \item Time points (time).  \item Unique patient
#' identifier (id).  \item FEV1-outcomes (y).  }
#' @references Mason D.P., Rajeswaran J., Li L., Murthy S.C., Su J.W.,
#' Pettersson G.B., Blackstone E.H. Effect of changes in postoperative
#' spirometry on survival after lung transplantation. \emph{J. Thorac.
#' Cardiovasc. Surg.}, 144:197-203, 2012.
#' @keywords datasets
#' @examples
#' data(spirometry, package = "boostmtree")
NULL






