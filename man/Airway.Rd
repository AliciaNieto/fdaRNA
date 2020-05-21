\name{Airway}
\alias{Airway}
\docType{data}
\title{
Airway smooth muscle cells
}
\description{
Read counts per gene for airway smooth muscle cell lines RNA-Seq experiment. The dataset contains four cell lines in two conditions: control and treatment with dexamethasone. It consists of 29391 RNA's that have more than one nonzero values. About half of the RNA's of the original data had two or more samples with a zero value and were discarded. 
}
\usage{data("Airway")}
\format{
  A data frame with 29391 observations on the following 8 variables.
  \describe{
    \item{\code{SRR1039508}}{a numeric vector}
    \item{\code{SRR1039509}}{a numeric vector}
    \item{\code{SRR1039512}}{a numeric vector}
    \item{\code{SRR1039513}}{a numeric vector}
    \item{\code{SRR1039516}}{a numeric vector}
    \item{\code{SRR1039517}}{a numeric vector}
    \item{\code{SRR1039520}}{a numeric vector}
    \item{\code{SRR1039521}}{a numeric vector}
  }
}
\details{
Variables:

Col. odd:  	SRR1039508, SRR1039512, SRR1039516 and SRR1039520. Untreated data: control group.

Col. even:  SRR1039509, SRR1039513, SRR1039517, SRR1039521. Treated data.
}
\source{
Bioconductor. Package 'airway'.
}
\references{
Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky- Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr, Tantisira KG, Weiss ST, Lu Q. ’RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.’ PLoS One. 2014 Jun 13;9(6):e99625. PMID: 24926665. GEO: GSE52778.

Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\examples{
outlier.detection(Airway, clases = 2, Cla = list(c(1,3,5,7),c(2, 4, 6 ,8)), MS = TRUE)
}
\keyword{datasets}