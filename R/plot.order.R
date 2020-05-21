#' Plots ordered datasets
#'
#' In the process of normalization, statistical depth based or quantile, each element of the sample 
#' is ordered in ascending order. This function plots the sample elements in this process.
#' 
#' @param ar Cell files from Affy micrarray experiment: a matrix where the columns are the sample 
#' elements, RNA-seq or microarray, and the rows the variables (the genes).
#' @param x A vector that is subset of the natural numbers up to the number of variables. 
#' It represents the domain, of the X-axis, to use for the zoomed plot. The default is the whole domain.
#' @param tp The type of plot that should be drawn. The default is "l" for lines.
#' @param xl A vector with the limits of the X-axis. The default is the range of "x".
#' @param yl A vector with the limits of the Y-axis. The default is the range of "ar" when restricted 
#' to the domain "x".
#' @param xlb X-axis label for the  multidimensional scaling and/or spectral map plots. The default 
#' is empty. It is only used if "MS" or/and "SM" is/are set to TRUE.
#' @param ylb Y-axis label for the  multidimensional scaling and/or spectral map plots. The default 
#' is empty. It is only used if "MS" or/and "SM" is/are set to TRUE.
#' @param m Main title plot. The default is empty
#' @param ... Additional arguments to be passed to the multidimensional scaling and/or spectral map 
#' plots, such as graphical parameters. It is only used if "MS" or/and "SM" is/are set to TRUE.
#'
#' @return does not return
#' 
#' @details
#' This function plots allows to plot a zoom of the domain. The color shows the deepness of the 
#' sample elemnt with respect to the sample, using \emph{prof.funct}. The colours go from blue 
#' meaning low depth to red meaning high depth, through green and yellow.
#'
#' @export
#' @importFrom grDevices rainbow
#' @importFrom graphics matlines
#'
#' @author Alicia Nieto-Reyes and Javier Cabrera
#'
#' @keywords robust multivariate non-parametric.
#' 
#' @seealso \code{\link{normalization}} and \code{\link{prof.funct}}
#' 
#' @references 
#' Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene 
#' expression data. Preprint.
#'
#' @examples
#' # Plots the Sialin data set after 6 hours of gestation, in the ascending order used 
#' # in the process of normalization
#' plot_order(Sialin6)
#' 
#' # A zoom of previous plot in the interval [20000 , 30000]
#' plot_order(Sialin6, c(20000:30000))
#'
plot_order <- function(ar,x=1:nrow(ar),tp='l', xl=range(x), yl=c(min(ar[x,]), max(ar[x,])), xlb='',ylb='',m='',...){
  sar = apply(ar, 2, sort); O=ncol(sar); 
  r=grDevices::rainbow(O, s = 1, v = 1, start = 0, end = .68, alpha = 1)
  po=prof.funct(sar)[[1]]*O; 
  plot(x, sar[x,1],type=tp,xlim=xl,ylim=yl,xlab=xlb,ylab=ylb,main=m,...)
  for (j in 1:O){w=which(po==j); if(length(w)!=0) {graphics::matlines(x,sar[x,w], type='l',col=r[O-j+1])}}
}
