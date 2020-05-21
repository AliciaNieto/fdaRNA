\name{Tissue}
\alias{Tissue}
\docType{data}
\title{
Three treatments to tissue samples data frame 
}
\description{
The Tissue data set consists of the study of the tissue, on mice, using three different treatments. It has 3,487 genes and 41 microarrays. Treatment one and two contain each 11 microarrays each and treatment three 19 microarrays. This is raw data.
}
\usage{data("Tissue")}
\format{
  A data frame with 3487 observations on the following 41 variables.
  \describe{
    \item{\code{X1}}{a numeric vector}
    \item{\code{X2}}{a numeric vector}
    \item{\code{X3}}{a numeric vector}
    \item{\code{X4}}{a numeric vector}
    \item{\code{X5}}{a numeric vector}
    \item{\code{X6}}{a numeric vector}
    \item{\code{X7}}{a numeric vector}
    \item{\code{X8}}{a numeric vector}
    \item{\code{X9}}{a numeric vector}
    \item{\code{X10}}{a numeric vector}
    \item{\code{X11}}{a numeric vector}
    \item{\code{X12}}{a numeric vector}
    \item{\code{X13}}{a numeric vector}
    \item{\code{X14}}{a numeric vector}
    \item{\code{X15}}{a numeric vector}
    \item{\code{X16}}{a numeric vector}
    \item{\code{X17}}{a numeric vector}
    \item{\code{X18}}{a numeric vector}
    \item{\code{X19}}{a numeric vector}
    \item{\code{X20}}{a numeric vector}
    \item{\code{X21}}{a numeric vector}
    \item{\code{X22}}{a numeric vector}
    \item{\code{X23}}{a numeric vector}
    \item{\code{X24}}{a numeric vector}
    \item{\code{X25}}{a numeric vector}
    \item{\code{X26}}{a numeric vector}
    \item{\code{X27}}{a numeric vector}
    \item{\code{X28}}{a numeric vector}
    \item{\code{X29}}{a numeric vector}
    \item{\code{X30}}{a numeric vector}
    \item{\code{X31}}{a numeric vector}
    \item{\code{X32}}{a numeric vector}
    \item{\code{X33}}{a numeric vector}
    \item{\code{X34}}{a numeric vector}
    \item{\code{X35}}{a numeric vector}
    \item{\code{X36}}{a numeric vector}
    \item{\code{X37}}{a numeric vector}
    \item{\code{X38}}{a numeric vector}
    \item{\code{X39}}{a numeric vector}
    \item{\code{X40}}{a numeric vector}
    \item{\code{X41}}{a numeric vector}
  }
}
\details{
Variables:

Col 1-11:   treat1.X1 through treat1.X11. Treatment 1 

Col 12-22:  treat2.X1 through treat2.X11. Treatment 2 

Col 23-41:  treat3.X1 through treat3.X19. Treatment 3 
}
\source{
Own
}
\references{
Amaratunga D, Cabrera J, Shkedy Z. Exploration and analysis of DNA microarray and other high dimensional data. J. Wiley & Sons, 2014. 

Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\examples{
normalization.boxplot(Tissue, before = TRUE, par.r = TRUE)
}
\keyword{datasets}
