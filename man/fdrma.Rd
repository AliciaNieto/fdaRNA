\name{fdrma}
\alias{fdrma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional Data Robust Microarray Analysis
}
\description{
This function performs robust microarray analysis for
a sample of microarray cell files from a microarray
experiment with Affymetric chips.
}
\usage{
fdrma(slc0cel, method = c("fdn", "myrma"), mod = c("MP", "MM", "LS"), no.cores = floor(detectCores()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slc0cel}{
Cell files from Affy micrarray experiment}
  \item{method}{ Normalization method: either "fdn" (funcional data normalization) or "qn" quantile normalization.

}
  \item{mod}{ one of the following:
 "MP"= median polish, "MM"= Tukey's biweight M-estimator,"M"= Huber M-estimator,"LS"= Least squares.
}
  \item{no.cores}{
default is floor(detectCores())-1
}
}
\details{
This function perform robust microarray analysis with the options of using
functional data  normalization instead of quantile  normalization.
In addition the user may specify a robust regression M-estimator of regression
("M" or "MM") instead of one step median polish that was used in the original RMA.

MM-estimator is an M-estimator that uses the psi function given by Tukey's biweight function. The M-estimator uses Huber's psi function.
}
\value{
This function returns an array of gene expressions with G genes by n samples, that summarizes the probes of the original datasets.}
\references{
Nieto A., Cabrera J. (2020) Statistical depth Based Normalization and Outlier Detection of Gene Expression Data.(to appear)

Amaratunga, D.; Cabrera, J. (2001). "Analysis of Data from Viral DNA Microchips. Journal of the American Statistical Association. 96 (456): 1161.

Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data. Nucleic Acids Research 31(4)

}
\author{
A. Nieto and J. Cabrera
}
\examples{
library(affy)
# Path to file with cel files from experiment
fns <- list.celfiles(path="cel_filesday0",full.names=TRUE)

# Read Affy files into R
slc0cel <- ReadAffy(filenames=fns)

# Get rma file and fdn files with gene expressions
slc0.fdn= fdrma(slc0cel, method = "fdn", mod = "MM")
slc0.rma= fdrma(slc0cel, method = "rma")

# Compare the functional forms of the gene expressions
plot(sort(slc0.rma[,1]),sort(slc0.fdn[,1]),type="l",col=2)
}

\keyword{ robust }
\keyword{ multivariate }
