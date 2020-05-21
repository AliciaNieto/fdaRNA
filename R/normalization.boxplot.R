#' Boxplot of gene expressions
#'
#' Boxplot of the "log(ar+1)" transform of the data before and after the statistical functional 
#' depth based normalization. This is commonly carried out to show the adequacy of the procedure.
#'  
#' @param ar Cell files from Affy micrarray experiment: a matrix where the columns are the sample 
#' elements, RNA-seq or microarray, and the rows the variables (the genes).
#' @param namesN If TRUE, the labels of each RNA-seq's or microarray are substituted by a number 
#' corresponding to its position in data matrix. The default is TRUE.
#' @param before If TRUE, the "log(ar+1)" transform of the data before the normalization step is
#' plotted. The default is FALSE.
#' @param after If TRUE, the "log(ar+1)" transform of the data after the normalization step is plotted. 
#' The default is TRUE.
#' @param par.r If TRUE, and both "before" and "after" are also TRUE, the "before" and "after" boxplots 
#' are displayed in a row. The default is FALSE.
#' @param par.c IIf TRUE, and both "before" and "after" are also TRUE, the "before" and "after" boxplots 
#' are displayed in a column. The default is FALSE.
#' @param arg1 Overall title for the boxplot before the normalization step. The default is \emph{Before}.
#' @param arg2 Overall title for the boxplot after the normalization step. The default is \emph{After}.
#'
#' @return The boxplot of the "log( +1)" transform of the data before and/or after the normalization step.
#'
#' @details
#' Boxplot of the "log(ar+1)" transform of the data before and after the normalization step. This is 
#' commonly carried out to show the adequacy of the procedure.
#'
#' @export
#' @importFrom graphics par
#'
#' @author A. Nieto and J. Cabrera
#'
#' @keywords robust multivariate non-parametric.
#' 
#' @references 
#' Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression 
#' data. Preprint.
#'
#' @examples
#' ## Computes the normalization of the RNA-seq's in the Airway dataset.
#' N = normalization(Airway)
#' 
#' ## Displays in a row the boxplots of the "log(+1)" transform 
#' ## of the data before and after normalization.
#' 
#' normalization.boxplot(Airway, before = TRUE, par.r = TRUE)
#'
normalization.boxplot <-
function(ar, namesN=T, before=F, after=T, par.r=F, par.c=F, arg1="Before", arg2="After"){
  if(namesN){names(ar)<-1:dim(ar)[2]}
  if(before){
    if(after){
      if(par.r){graphics::par(mfrow=c(1,2))}
      if(par.c){graphics::par(mfrow=c(2,1))}
    }
    boxplot(log(ar+1),main=arg1)}
  if(after){boxplot(log(normalization(ar)+1),main=arg2)}
  par(mfrow=c(1,1))
}
