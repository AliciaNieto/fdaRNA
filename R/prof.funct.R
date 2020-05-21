#' Emprirical statistical depth
#'
#' Computes the sample statistical depth of a given sample using the procedure described in 
#' Nieto-Reyes (2011) and Nieto-Reyes and Cabrera (2020).
#' 
#' @param x A matrix where the columns are the sample elements and the rows the variables. 
#' The matrix can represent functional or multivariate data. It requires at least two variables. 
#' When applied to gene expression data, each column is a RNA-seq or microarray and the rows 
#' represent the genes. 
#'
#' @return The function returns a list containing the following components:
#' \itemize{
#'    \item{depth.va}{A vector with the depth values of the sample elements.}
#'    \item{deepest.ele}{The identification number (column) of the deepest(s) elements.}
#'    \item{distant.val}{In the process of computing the depth values, these are ordered in pairs. 
#'    There is a number of pairs equal to the integer part of the sample size divided by two. This 
#'    vector contains the distance among each two elements of a pair, in decreasing order.}
#' }
#' 
#' @export
#'
#' @author Alicia Nieto-Reyes and Javier Cabrera
#'
#' @keywords robust multivariate non-parametric.
#' 
#' @seealso \code{\link{normalization}} and \code{\link{prof.funct}}
#' 
#' @references 
#' Nieto-Reyes A. (2011) On the Properties of Functional Depth. In: Ferraty F. (eds) 
#' Recent Advances in Functional Data Analysis and Related Topics. Contributions to Statistics. 
#' Physica-Verlag HD.
#' 
#' Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of 
#' gene expression data. Preprint.
#'
#' @examples
#' # Applies the "prof.funct" function to the Tissue dataset
#' p = prof.funct(Tissue)
#' 
#' # Gives the depth values of each of the 41 microarrays in the Tissue dataset
#' p$depth.val
#' 
#' # Gives the deepest microarray of the Tissue dataset, the one with highest depth value.
#' p$deepest.ele
#'
prof.funct <-function(x){
  mg=as.matrix(dist(t(x))) 
  m=mg; k=ncol(x); for (i in 1:k) {m[i,i:k]=0}
  r=matrix(0,k,k); c=matrix(0,k,k); for (i in 1:k){r[i,]=i; c[,i]=i}
  p=rep(0,k); e=0
  if((k %% 2)==0) {fu=k-2; di=rep(0,k/2)} else {fu=k-1; di=rep(0,fu/2)} 
  cont=0
  for (i in which(((1:fu) %% 2)!=0)) 
  {cont=cont+1; di[cont]=max(m)
    w=which.max(m); ro=r[w][1]; co=c[w][1]; u=min(mg[-c(e,ro),ro]); du=min(mg[-c(e,co),co])
    if (u>du) 
    {p[ro]=i/k; p[co]=(i+1)/k} 
    else 
    {
      if (u==du) 
      {s=sample(c(1,0)); p[ro]=(i+s[1])/k; p[co]=(i+s[2])/k} 
      else 
      {p[ro]=(i+1)/k; p[co]=i/k}
    }
    m[ro,]=0; m[,co]=0; m[,ro]=0; m[co,]=0; e=c(e,ro,co)
  }
  p[which(p==0)]=1; if((k %% 2)==0){di[cont+1]=max(m)}
  list(depth.val = p, deepest.ele = which(p==max(p)), distant.val = di)
}
