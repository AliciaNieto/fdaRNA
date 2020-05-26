#' Functional Data Robust Microarray Analysis
#'
#' This function performs robust microarray analysis for
#' a sample of microarray cell files from a microarray
#' experiment with Affymetric chips.
#'
#' @param slc0cel Cell files from Affy micrarray experiment.
#' @param method  Normalization method: either "fdn" (funcional data normalization) or 
#' "qn" quantile normalization.
#' @param mod one of the following: "MP"= median polish, 
#' "MM"= Tukey's biweight M-estimator,"M"= Huber M-estimator,"LS"= Least squares.
#' @param no.cores default is floor(detectCores())-1
#'
#' @return This function returns an array of gene expressions with G genes by n samples, 
#' that summarizes the probes of the original datasets.
#'
#' @details
#' functional data  normalization instead of quantile  normalization. In addition the user may 
#' specify a robust regression M-estimator of regression ("M" or "MM") instead of one step median 
#' polish that was used in the original RMA. MM-estimator is an M-estimator that uses the psi function 
#' given by Tukey's biweight function. The M-estimator uses Huber's psi function.
#'
#' @export
#' @import parallel MASS stats
#' @importFrom affy pm bg.correct.rma
#' @importFrom doSNOW registerDoSNOW
#'
#' @author A. Nieto and J. Cabrera
#'
#' @keywords robust multivariate
#'
#' @examples
#' \dontrun{
#' library(affy)
#' # Path to file with cel files from experiment
#' fns <- list.celfiles(path="cel_filesday0",full.names=TRUE)
#' # Read Affy files into R
#' slc0cel <- ReadAffy(filenames=fns)
#' 
#' # Get rma file and fdn files with gene expressions
#' slc0.fdn= fdrma(slc0cel, method = "fdn", mod = "MM")
#' slc0.rma= fdrma(slc0cel, method = "rma")
#' # Compare the functional forms of the gene expressions
#' plot(sort(slc0.rma[,1]),sort(slc0.fdn[,1]),type="l",col=2)
#' }
#' 
fdrma <-function(slc0cel,method=c("fdn","myrma"), mod=c("MP","MM","LS"),no.cores = floor(detectCores())) {
   # require(affy); require(MASS)
  if(method[1]=="myrma") { my.PM=affy::pm(normalize.AffyBatch.quantiles(slc0cel)) } else
    if(method[1]=="fdn"){
      slc0celbg=affy::bg.correct.rma(slc0cel)
      pmslc0 = affy::pm(slc0celbg)
      #    pmslc0 <- apply(pmslc0, 2, bg.adjust)
      my.PM= normalization(pmslc0)
      z= floor (min((my.PM))); if(z<0) my.PM=my.PM-z
      #    for(i in 1:ncol(my.PM)) my.PM[,i]= my.PM[rank(pmslc0[,i]),i]
    }
  pG = affy::probeNames(slc0cel)
  pn = unique(pG)
  
  G = length(pn)
  n = ncol(my.PM)
  xx = log2(my.PM)
  
  cl = parallel::makeCluster(no.cores)
  doSNOW::registerDoSNOW(cl)
  z = split(data.frame(xx),pG)
  fun1=function(x){ xxx = medpolish (as.matrix(x),trace.iter=F);xxx$overall+xxx$col }

       fun3 = function(x) {
      i0 = nrow(x); nc=ncol(x)
      tt11 = data.frame( y=c(as.matrix(x)), x=factor(rep(colnames(x),rep(i0,nc) )),
                         p=factor(rep(1:i0,nc)))
      tt3 = lm(y~p+x,data=tt11)$coef
      tt3[1]+c(0,tt3[-(1:i0)])+mean(c(0,tt3[2:i0]))
    }
   
  fun2=function(x) {
    i0 = nrow(x); nc=ncol(x)
    tt11 = data.frame( y=c(as.matrix(x)), x=factor(rep(colnames(x),rep(i0,nc) )),
                       p=factor(rep(1:i0,nc)))
    for(i1 in 1:1000) {
      tt4=try( MASS::rlm(y~p+x,data=tt11,method="MM",init="ls")$coef)
      if(is.numeric(tt4[1])) break
    }
    
    tt4[1]+c(0,tt4[-(1:i0)])+median(c(0,tt4[2:i0]))
  }
  fun= if(mod=="MP") fun1 else if(mod=="MM") fun2 else fun3
  zz3=parallel::parSapply(cl,z,fun)
  parallel::stopCluster(cl)
  t(zz3)
}
