#' Statistical depth based normalization
#'
#' Gene expression normalization equates the scales at which the various gene expressions 
#' have been measured. This function follows the procedure in Nieto-Reyes and Cabrera (2020) 
#' that equates to the statistical sample functional depth in Nieto-Reyes (2011).
#' 
#' @param x Cell files from Affy micrarray experiment: a matrix where the columns are the sample 
#' elements, RNA-seq or microarray, and the rows the variables (the genes).
#'
#' @return It returns a matrix of the same dimension as the argument, which consists of the
#' normalized data.
#'
#' @export
#'
#' @author A. Nieto and J. Cabrera
#'
#' @keywords robust multivariate
#' 
#' @references 
#' Amaratunga D, Cabrera J, Shkedy Z. Exploration and analysis of DNA microarray and other high 
#' dimensional data. J. Wiley & Sons, 2014. 
#' Nieto-Reyes A. (2011) On the Properties of Functional Depth. In: Ferraty F. (eds) 
#' Recent Advances in Functional Data Analysis and Related Topics. Contributions to Statistics. 
#' Physica-Verlag HD.
#' Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier 
#' detection of gene expression data. Preprint.
#'
#' @examples
#' ## Computes the normalization of the RNA-seq's in the Airway dataset.
#' N = normalization(Airway)
#' 
#' ## Displays in a row the boxplots of the "log( +1)" transform of 
#' ## the data before and after normalization.
#' normalization.boxplot(Airway, before = TRUE, par.r = TRUE)
#' 
normalization <-function (x){
   sar = apply(x, 2, sort); pom=prof.funct(sar)[[2]];
  if (length(pom)>1) xxm=sort(rowMeans(x[,pom])) else xxm=sort(x[,pom])
  xr <- c(apply(x, 2, rank))
  array(approx(1:nrow(x), xxm, xr)$y, dim(x), dimnames(x))}
