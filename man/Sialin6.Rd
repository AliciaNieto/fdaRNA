\name{Sialin6}
\alias{Sialin6}
\docType{data}
\title{
Embryo stage: S1c17A5 gene knocked out versus wild type
}
\description{
Compares the gene expression profiles of 6 knocked out mice, mice for which the S1c17A5 gene has been knocked out, versus 6 wild type mice. That is, 12 microarrays divided in two groups of six, with 45,101 genes each. This is a microarray experiment that was repeated at two time points, one in day first of development during the embryo stage (after six hours).
The other, that can be found under Sialin18, at 18 days of age during the adult stage of the mice. he  biology suggests that there should be small differences after just a few hours of gestation and larger differences after 18 days.
}
\usage{data("Sialin6")}
\format{
  A data frame with 45101 observations on the following 12 variables.
  \describe{
    \item{\code{X2760.CEL}}{a numeric vector}
    \item{\code{X2761.CEL}}{a numeric vector}
    \item{\code{X2762.CEL}}{a numeric vector}
    \item{\code{X2763.CEL}}{a numeric vector}
    \item{\code{X2764.CEL}}{a numeric vector}
    \item{\code{X2765.CEL}}{a numeric vector}
    \item{\code{X2766.CEL}}{a numeric vector}
    \item{\code{X2767.CEL}}{a numeric vector}
    \item{\code{X2768.CEL}}{a numeric vector}
    \item{\code{X2769.CEL}}{a numeric vector}
    \item{\code{X2770.CEL}}{a numeric vector}
    \item{\code{X2771.CEL}}{a numeric vector}
  }
}
\details{
Variables:
  
Col. 1,4,6,7,9,10:    X2760.CEL, X2763.CEL, X2765.CEL, X2766.CEL, X2768.CEL, X2769.CEL. Group 1

Col. 2,3,5,8,11,12:   X2761.CEL, X2762.CEL, X2764.CEL, X2767.CEL, X2770.CEL, X2771.CEL. Group 2
}
\source{
Own
}
\references{
Amaratunga D, Cabrera J, Shkedy Z. Exploration and analysis of DNA microarray and other high dimensional data. J. Wiley & Sons, 2014. 

Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\examples{
plot_order(Sialin6, m='Sialin after 6 hours dataset')
}
\keyword{datasets}
