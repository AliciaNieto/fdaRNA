\name{Sialin18}
\alias{Sialin18}
\docType{data}
\title{
Adult stage: S1c17A5 gene knocked out versus wild type
}
\description{
Compares the gene expression profiles of 6 knocked out mice, mice for which the S1c17A5 gene has been knocked out, versus 6 wild type mice. That is, 12 microarrays divided in two groups of six, with 45,101 genes each. This is a microarray experiment that was repeated at two time points, one in day first of development during the embryo stage, that can be found under Sialin6.
This one, at 18 days of age during the adult stage of the mice. he  biology suggests that there should be small differences after just a few hours of gestation and larger differences after 18 days.
}
\usage{data("Sialin18")}
\format{
  A data frame with 45101 observations on the following 12 variables.
  \describe{
    \item{\code{X2720.CEL}}{a numeric vector}
    \item{\code{X2721.CEL}}{a numeric vector}
    \item{\code{X2722.CEL}}{a numeric vector}
    \item{\code{X2723.CEL}}{a numeric vector}
    \item{\code{X2724.CEL}}{a numeric vector}
    \item{\code{X2725.CEL}}{a numeric vector}
    \item{\code{X2726.CEL}}{a numeric vector}
    \item{\code{X2727.CEL}}{a numeric vector}
    \item{\code{X2728.CEL}}{a numeric vector}
    \item{\code{X2729.CEL}}{a numeric vector}
    \item{\code{X2730.CEL}}{a numeric vector}
    \item{\code{X2731.CEL}}{a numeric vector}
  }
}
\details{
Variables:
  
Col. 1,4,6,7,9,10:    X2720.CEL, X2723.CEL, X2725.CEL, X2726.CEL, X2728.CEL, X2729.CEL. Group 1

Col. 2,3,5,8,11,12:   X2721.CEL, X2722.CEL, X2724.CEL, X2727.CEL, X2730.CEL, X2731.CEL. Group 2   
}
\source{
Own
}
\references{
Amaratunga D, Cabrera J, Shkedy Z. Exploration and analysis of DNA microarray and other high dimensional data. J. Wiley & Sons, 2014. 

Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\examples{
plot_order(Sialin18, x=20000:21000, m='Zoom: Sialin after 18 days dataset')
}
\keyword{datasets}
