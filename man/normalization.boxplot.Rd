\name{normalization.boxplot}
\alias{normalization.boxplot}
\title{
Boxplot of gene expressions
}
\description{
Boxplot of the "log(ar+1)" transform of the data before and after the statistical functional depth based normalization. This is commonly carried out to show the adequacy of the procedure.
}
\usage{
normalization.boxplot(ar, namesN = T, before = F, after = T, par.r = F, par.c = F, arg1 = "Before", arg2 = "After")
}
\arguments{
  \item{ar}{
  Cell files from Affy micrarray experiment: a matrix where the columns are the sample elements, RNA-seq or microarray, and the rows the variables (the genes).
}
  \item{namesN}{
If TRUE, the labels of each RNA-seq's or microarray are substituted by a number corresponding to its position in data matrix. The default is TRUE.
}
  \item{before}{
If TRUE, the "log(ar+1)" transform of the data before the normalization step is plotted. The default is FALSE.
}
  \item{after}{
If TRUE, the "log(ar+1)" transform of the data after the normalization step is plotted. The default is TRUE.
}
  \item{par.r}{
If TRUE, and both "before" and "after" are also TRUE, the "before" and "after" boxplots are displayed in a row. The default is FALSE.
}
  \item{par.c}{
If TRUE, and both "before" and "after" are also TRUE, the "before" and "after" boxplots are displayed in a column. The default is FALSE.
}
  \item{arg1}{
Overall title for the boxplot before the normalization step. The default is $Before$.
}
  \item{arg2}{
Overall title for the boxplot after the normalization step. The default is $After$.
}
}
\details{
Boxplot of the "log(ar+1)" transform of the data before and after the normalization step. This is commonly carried out to show the adequacy of the procedure.
}
\value{
The boxplot of the "log( +1)" transform of the data before and/or after the normalization step.
}
\references{
Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\author{
Alicia Nieto-Reyes and Javier Cabrera
}
\seealso{
\code{\link{normalization}}
}
\examples{
## Computes the normalization of the RNA-seq's in the Airway dataset.
N = normalization(Airway)

## Displays in a row the boxplots of the "log( +1)" transform of the data before and after normalization.
normalization.boxplot(Airway, before = T, par.r = T)

}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
